#include <cstdint>
#include <vector>
#include <iostream>
#include <algorithm>

#include <eigen3/Eigen/Dense>       // cross, dot product.

#include "MathTypes.hpp"
#include "LevelSet.hpp"
#include "InitialShapes.hpp"
#include "common/Utils.hpp"

using namespace FEM;


namespace {
    // Level set function reset weights (see update()) for phi_ and phi_ifn_.
    static const std::pair<double, double> PHI_RESET_WEIGHTS(0.75, 0.25); 
}



int LevelSet::init( const Parameters& par )
{
    par_ = par;

    std::cout << "  - Creating mesh..." << std::flush;
    if (mesh_.createRegularMesh( par_, nThreads_ )) {
        std::cerr << "Error: Failed to create a mesh." << std::endl;
        return -1;
    }
    std::cout << "done." << std::endl;


    std::cout << "  - Assembling advection mass... " << std::flush;
    uint32_t n = mesh_.p().size();
    M_ = sparse_matrix(n,n);
    auto wM = WEAK_FORM( u.cwiseProduct(v) )
    assembler_.assemble_bilin( mesh_, wM, {}, M_ );
    std::cout << "done." << std::endl;


    std::cout << "  - Setting initial domain shape..." << std::flush;
    phi_ = VecXd::Zero(n);
    phi_ifn_ = phi_;
    Vec3f dim( {par_["Axis_a"], par_["Axis_b"], par_["Axis_c"]} );

    switch (int s = (int)par_["Initial shape"]) {
    case ELLIPSOID:
        Ellipsoid( mesh_, dim[0], dim[1], dim[2], phi_ );
        break;
    case CUBOID:
        Cuboid( mesh_, dim[0], dim[1], dim[2], phi_ );
        break;
    case LEVEQUE:
        LeVeque( mesh_, phi_ );
        break;
    default:
        std::cerr << "Unrecognized value of 'Initial shape': " << s << std::endl;
        return -1;
    }
    std::cout << "done." << std::endl;


    // Upload mesh nodes, initial phi to a CUDA device.
    // Stores the pointer only if no GPU available.
    std::cout << "  - Uploading vertices to GPU..." << std::flush;
    dist_.init( mesh_.p(0).data(), phi_.data(), n, nThreads_ );
    std::cout << "done." << std::endl;

    // Compute initial distance function. To given smoothes possible initial shapes,
    // don't preserve the pre-reset values
    resetPhi_();
    phi_ = phi_ifn_;

    // Initialize the solver.
    double CG_tol = par_["Advection tolerance"];
    if (CG_tol > 0.0)
        solver_.setTolerance( CG_tol );
    int CG_iter = (int)par_["Max. CG iterations"];
    if (CG_iter > 0)
        solver_.setMaxIterations( CG_iter );

    return 0;
}



int LevelSet::update( const VecXd& flow )
{
    double dt = par_["Advection dt"];
    if (dt < FLT_MIN || flow.norm() < FLT_MIN)
        return 0;
 
    //
    // Solve advection on phi (updates interface position).
    //
    uint32_t n = mesh_.p().size();
    if (flow.size() != 4*n) {
        std::cerr << __FUNCTION__ << "(): Stokes flow vector doesn't match "
                  << " with the mesh." << std::endl;
        return -1;
    }

    chrono_time timeStart;
    tic( timeStart );
    sparse_matrix A(n,n);
    auto wA = WEAK_FORM( V1.cwiseProduct(ux.cwiseProduct(v)) +
                         V2.cwiseProduct(uy.cwiseProduct(v)) +
                         V3.cwiseProduct(uz.cwiseProduct(v)) )
    auto& V = flow.head(3*n);         // 3D velocity, ignore pressure
    if (V.hasNaN()) {
        std::cerr << __FUNCTION__ << "(): Take your NaNsense elsewhere!" << std::endl;
        return -1;
    }

    assembler_.assemble_bilin( mesh_, wA, V, A );
    std::cout << "  advection assembly took " <<  toc( timeStart ) << " ms." 
              << std::endl;
    

    tic( timeStart );
    // Implicit Euler
    // solver_.compute( M_ + dt*A );
    // VecXd phi_n = solver_.solveWithGuess( M_ * phi_, phi_ );
    // Crank-Nicolson
    solver_.compute( M_ + 0.5*dt*A );
    VecXd phi_n = solver_.solveWithGuess( (M_ - 0.5*dt*A) * phi_, phi_ );

    std::cout << "  solver finished after " << toc( timeStart ) << " ms. and " 
              << solver_.iterations() << " iterations." << std::endl;
    std::cout << "  Error estimate: " << solver_.error() << std::endl;
    std::cout << "  phi delta norm " << (phi_-phi_n).norm() << std::endl;

    if (solver_.info() != Eigen::Success) {
        std::cerr << __FUNCTION__ << "(): Solving advection equation failed "
                  << "for unknown reason." << std::endl;
        return -1;
    }

    phi_.swap( phi_n );

    //
    // Recompute phi to make sure it is true distance function.
    //
    std::cout << "  Reseting level set... " << std::endl;
    resetPhi_();
    phi_ = PHI_RESET_WEIGHTS.first*phi_ + PHI_RESET_WEIGHTS.second*phi_ifn_;

    return 0;
}



int LevelSet::interpolatePhi( const Vec3f& p, float& phi, Vec3f& normal ) const
{
    Vec4f bc;   // barycentric coordinates of p inside t
    Vec4i t;    // p enclosing tetrahedron
    if (mesh_.getEnclosingTetrahedron(p, bc, t)) {
        std::cerr << __FUNCTION__ << ": Failed to find element containing "
                  << "knot " << p.transpose() << std::endl;
        return -1;
    }

    // interpolate level set normal at p
    normal = Vec3f( {0.0f, 0.0f, 0.0f} );
    for (int j=0; j<4; j++)
        normal += bc(j) * normals_[t(j)] / normals_[t(j)].norm();
    normal /= normal.norm();

    // interpolate phi at p
    Vec4d phi_t( {phi_(t(0)), phi_(t(1)), phi_(t(2)), phi_(t(3))} );
    phi = bc.dot( phi_t.cast<float>() );

    return 0;
}



std::vector<int> LevelSet::getInteriorElements() const
{
#pragma omp declare reduction \
        (merge : std::vector<int> : omp_out.insert(omp_out.end(), \
                                    omp_in.begin(), omp_in.end()))
    std::vector<int> elements;
#pragma omp parallel for num_threads( nThreads_ ) reduction( merge: elements )

    for (size_t i=0; i<mesh_.t().size(); i++) {
        auto& t = mesh_.t(i);
        if ( phi_(t(0)) < 0.0 && phi_(t(1)) < 0.0 && phi_(t(2)) < 0.0 &&
             phi_(t(3)) < 0.0 )
            elements.push_back(i);
    }

    return elements;
}



std::vector<int> LevelSet::getInteriorNodes() const
{
    std::vector<int> nodes;
    nodes.reserve(phi_.size());

    for (size_t i=0; i<mesh_.p().size(); i++) {
        if (phi_(i) < 0.0)
            nodes.push_back(i);
    }

    return nodes;
}



std::vector<float> LevelSet::getInterfaceCurvature() const
{
    //
    // Curvature at mesh nodes.
    //
    std::vector<Vec3f> g;
    assembler_.grad_nodes( mesh_, phi_ifn_, g );

    size_t n = g.size();
    VecXd g_x(n), g_y(n), g_z(n);
    for (size_t i=0; i<n; i++) {
        Vec3f gn = g[i] / g[i].norm();
        g_x[i] = gn(0);                   // (dphi / dx) / |nabla phi|
        g_y[i] = gn(1);                   // (dphi / dy) / |nabla phi|
        g_z[i] = gn(2);
    }

    std::vector<Vec3f> gg_x, gg_y, gg_z;
    assembler_.grad_nodes( mesh_, g_x, gg_x );  // nabla (dphi / dx) / |nabla phi|
    assembler_.grad_nodes( mesh_, g_y, gg_y );
    assembler_.grad_nodes( mesh_, g_z, gg_z );

    std::vector<float> k_mesh(n);
    for (size_t i=0; i<n; i++)
        k_mesh[i] = gg_x[i](0) + gg_y[i](1) + gg_z[i](2);   // divergence

    //
    // Curvature interpolated at the interface nodes along mesh edges.
    //
    std::vector<float> k_if;
    k_if.reserve( e2ifp_.size() ); 
    for (auto& m : e2ifp_) {
        Vec2i e = mesh_.edges( m.first );
        float t = phi_[e(0)] / (phi_[e(0)] - phi_[e(1)]);
        k_if.push_back( (1.0-t)*k_mesh[e(0)] + t*k_mesh[e(1)] );
    }

    return k_if;
}


//
// Private methods
//


void LevelSet::findInterfaceNodes_()
{
    ifp_.clear();
    ifp_.reserve( mesh_.edges().size() );
    ifp2e_.clear();
    ifp2e_.reserve( mesh_.edges().size() );
    e2ifp_.clear();         // std::map, hence no reserve
   
    //
    // Get edges crossing zero level set, construct list of interface points.
    //
    chrono_time timeStart;
    std::cout << "  - Finding interface-mesh intersection points..."; 
    tic( timeStart );

    Vec3d cm(0.0, 0.0, 0.0);        // center of (interface) mass

    for (size_t i=0; i<mesh_.edges().size(); i++)
    {
        auto& e = mesh_.edges(i);
        if (std::signbit( phi_(e[0]) ) == std::signbit( phi_(e[1]) ))
            continue;           // edge not crossing interface

        // If the edge nodes level set delta close to zero, then the interface
        // approximates the edge; setting the crossing point at halfway between
        // the edge nodes.
        double t = 0.5;
        double div = phi_[e(0)] - phi_[e(1)];
        if (std::abs(div) > FLT_MIN)
            t = phi_[e(0)] / div;

        Vec3d p0 = mesh_.p(e(0)).cast<double>();
        Vec3d p1 = mesh_.p(e(1)).cast<double>();
        Vec3d pif = (1.0-t)*p0 + t*p1;
        cm += pif;

        ifp_.push_back( pif.cast<float>() );
        ifp2e_.push_back(i);
        e2ifp_[i] = e2ifp_.size();
    }

    ifp_.shrink_to_fit();
    ifp2e_.shrink_to_fit();

    //
    // Compute interface mean norm and standard deviation from the center of mass.
    //
    cm /= ifp_.size();
    std::vector<double> ifNorm( ifp_.size() );
    for (size_t i=0; i<ifp_.size(); i++)
        ifNorm[i] = (ifp_[i].cast<double>() - cm).norm();

    double mean = std::accumulate(ifNorm.begin(), ifNorm.end(), 0.0) / ifNorm.size();
    double var = 0.0;
    for (auto v : ifNorm)
        var += (v-mean) * (v-mean);
    var /= ifNorm.size();

    std::cout << " done after " << toc( timeStart ) << " ms." << std::endl;
    std::cout << "    Interface norm mean " << mean << ", std. " << std::sqrt(var) 
              << std::endl;    
}



int LevelSet::reconstructInterface_()
{
    // Adds a new triangle with nodes n to ift (and updates p2t) with correct orientation; 
    // accumulates the triangle areas ifpArea_ for each node.
    auto addTriangle = [&]( const Vec3i& n, int i, auto& ift, auto& p2t )
    {
        Vec3f e1 = ifMesh_.p(n(1)) - ifMesh_.p(n(0));
        Vec3f e2 = ifMesh_.p(n(2)) - ifMesh_.p(n(0));
        Vec3f normal = e1.cross(e2);

        // Find the reference node that doesn't coincide with the interface for 
        // triangle orienting.
        int p_ref = mesh_.t(i)(0);
        int k = 0;
        while (std::fabs(phi_(p_ref)) <= FLT_MIN)
            p_ref = mesh_.t(i)(++k);

        Vec3f e_ref = mesh_.p(p_ref) - ifMesh_.p(n(0));
        if (std::signbit( normal.dot(e_ref) ) != std::signbit( phi_[p_ref] ))
            ift.push_back( {n(0), n(2), n(1), -1} );
        else
            ift.push_back( {n(0), n(1), n(2), -1} );

        for (int k=0; k<3; k++) {
            // Append the new triangle to the node -> triangles map.
            VecXi& ti = p2t[n(k)];
            VecXi ti_new(ti.size() + 1);
            ti_new << ti, ift.size()-1;
            p2t[n(k)] = ti_new;
            // Accumulate triangle areas.
            ifpArea_[n(k)] += 0.5 * normal.norm();
        }
    };

    chrono_time timeStart;
    std::cout << "  - Reconstructing interface..."; 
    tic( timeStart );

    ifMesh_.set_p(ifp_);            // add interface vertices to surface mesh
    size_t np = ifp_.size();
    ifpArea_ = std::vector<float>(np, 0.0f);
    std::vector<Vec4i, Eigen::aligned_allocator<Vec4i>> ift;    // surface triangles
    ift.reserve(2*np);                                           // 2*np a rough estimate
    std::vector<VecXi> p2t(np, VecXi(0));

    ifElements_.clear();
    ifElements_.reserve(2*np);

    // For elements at the interface there's either 3 or 4 edges crossing
    // the interface. For 3 edges a single triangle is added; for 4 edges two
    // triangles.
    for (size_t i=0; i<mesh_.t().size(); i++)
    {
        auto& t = mesh_.t(i);
        int s = std::signbit( phi_(t(0)) ) + std::signbit( phi_(t(1)) ) +
                std::signbit( phi_(t(2)) ) + std::signbit( phi_(t(3)) );  
        if (s == 4 || s == 0)
            continue;                   // t doesn't cross the interface

        std::vector<int> eIf;           // edges that cross the interface
        auto& t2e = mesh_.t2e(i);
        for (int j=0; j<t2e.rows(); j++) {
            auto& e = mesh_.edges( t2e[j] );
            if (std::signbit( phi_(e[0]) ) != std::signbit( phi_(e[1]) ))
                eIf.push_back( t2e[j] );
        }

        if (eIf.size() == 3) {
            addTriangle( { e2ifp_[eIf[0]], e2ifp_[eIf[1]], e2ifp_[eIf[2]] }, i, ift, p2t );
        }
        else if (eIf.size() == 4) {
            addTriangle( { e2ifp_[eIf[0]], e2ifp_[eIf[2]], e2ifp_[eIf[3]] }, i, ift, p2t );
            addTriangle( { e2ifp_[eIf[0]], e2ifp_[eIf[1]], e2ifp_[eIf[3]] }, i, ift, p2t );
        }
        else {
            std::cerr << "Failed to reconstruct implicit surface: Element with "
                      << eIf.size() << " edge crossings (only 3, 4 accepted), "
                      << "t2e must be incorrectly filled!" << std::endl;
            return -1;
        }

        ifElements_.push_back(i);
    }

    ifElements_.shrink_to_fit();
    ift.shrink_to_fit();
    ifMesh_.set_t(ift);
    ifMesh_.set_p2t(p2t);

    std::cout << " done after " << toc( timeStart ) << " ms." << std::endl;

    float meanArea = std::accumulate(ifpArea_.begin(), ifpArea_.end(), 0.0) / ifpArea_.size();
    std::cout << "    Mean interface patch area: " << meanArea << std::endl;
    std::cout << "    " << ifElements_.size() << " interface elements, " << ifp_.size() 
              << " interface nodes, " << ifMesh_.t().size() << " interface triangles."
              << std::endl; 


    return 0;
}



void LevelSet::resetPhi_()
{
    findInterfaceNodes_();
    reconstructInterface_();

    chrono_time timeStart;
    std::cout << "  - Computing interface/mesh distances..."; 
    tic( timeStart );

    VecXd phi_n = phi_;
    dist_.resetPhi( ifp_[0].data(), ifp_.size(), phi_n.data() );
    phi_ifn_ = phi_n;       // fully reset values
    
    // Retain/don't update the values for the nodes of interface crossing edges:
    for (int& e : ifp2e_) {
        Vec2i n = mesh_.edges(e);
        phi_n(n(0)) = phi_(n(0));
        phi_n(n(1)) = phi_(n(1));
    }
    phi_ = phi_n;    
    std::cout << " done after " << toc( timeStart ) << " ms." << std::endl;

    std::cout << "  - Computing level set normals..."; 
    tic( timeStart );
    assembler_.grad_nodes( mesh_, phi_ifn_, normals_ );
    std::cout << " done after " << toc( timeStart ) << " ms." << std::endl;
}
