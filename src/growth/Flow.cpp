#include <iostream>
#include <cmath>

#include "Flow.hpp"
#include "ViscosityProfiles.hpp"
#include "common/Utils.hpp"

using namespace FEM;



namespace {
    // If true, uses the the previous solution as the initial guess for the 
    // iterative solver.
    const bool USE_INITIAL_GUESS = true;      
}



int Flow::init( const LevelSet& lset, const Parameters& par )
{   
    par_ = par;    
    uDeg_ = (int)par_["Stokes velocity degree"];
    if (uDeg_ < 1)
        return -1;

    auto& fem = lset.getAssembler();
    auto& mesh = lset.getMesh();
    size_t np = mesh.p().size();        // number of P1 dofs
    size_t nd = mesh.getNDofs(uDeg_);   // number of element dofs (P1 or P2)

    // 
    // Store the indices of nodes in which to solve the system.
    // 
    auto I = mesh.getInteriorDofIndices(uDeg_);
    size_t nI = I.size();
    sIdx_.resize(3*nI + np);
    for (size_t i=0; i<nI; i++) {
        sIdx_[i] = I[i];                // velocity x
        sIdx_[i+nI] = I[i] + nd;        // velocity y
        sIdx_[i+2*nI] = I[i] + 2*nd;    // velocity z
    }
    for (size_t i=0; i<np; i++)            
        sIdx_[i+3*nI] = i + 3*nd;       // pressure

    // Initialize the solution vector.
    up_ = VecXd(sIdx_.size());
    up_.setZero();

    //
    // Assemble the immutable matrices for pressure, velocity diverge (B*_)
    // pressure penalty (C_) and stabilization (D_) terms.
    //
    Bx_ = sparse_matrix(nd,np);
    auto wBx = WEAK_FORM( vx.cwiseProduct(u) )
    fem.assemble_bilin( mesh, wBx, {}, Bx_, uDeg_, 1 );

    By_ = sparse_matrix(nd,np);
    auto wBy = WEAK_FORM( vy.cwiseProduct(u) )
    fem.assemble_bilin( mesh, wBy, {}, By_, uDeg_, 1 );

    Bz_ = sparse_matrix(nd,np);
    auto wBz = WEAK_FORM( vz.cwiseProduct(u) )
    fem.assemble_bilin( mesh, wBz, {}, Bz_, uDeg_, 1 );

    M_ = sparse_matrix(np,np);
    auto wM = WEAK_FORM( u.cwiseProduct(v) )
    fem.assemble_bilin( mesh, wM, {}, M_, 1, 1 );

    N_ = sparse_matrix(nd,np);
    auto wN = WEAK_FORM( u.cwiseProduct(v) )
    fem.assemble_bilin( mesh, wN, {}, N_, uDeg_, 1 );

    //
    // Initialize the solver.
    //
    double CG_tol = par_["Stokes tolerance"];
    if (CG_tol > 0.0)
        solver_.setTolerance( CG_tol );
    int CG_iter = (int)par_["Max. CG iterations"];
    if (CG_iter > 0)
        solver_.setMaxIterations( CG_iter );    
    
    return 0;
}



int Flow::solve( const LevelSet& lset, const VecXd& source, VecXd& flow )
{
    if (!USE_INITIAL_GUESS)
        flow.setZero();
    return Stokes_( lset, source, flow );
}


//
//  PRIVATE METHODS
//


void Flow::constructMatrices_( const LevelSet& lset, const VecXd& mu )
{
    // If constant viscosity and A_ has already been assembled, nothing to do:
    if (par_["Viscosity interior"] - par_["Viscosity exterior"] <= FLT_MIN && 
        A_.nonZeros() > 0)
        return;
 
    auto& fem = lset.getAssembler();
    auto& mesh = lset.getMesh();
    size_t np = mesh.p().size();          // number of P1 dofs
    size_t nd = mesh.getNDofs(uDeg_);   // number of dofs (P1 or P2)

    // Frees up memory used by given sparse matrices.
    auto freeSparseMemory_ = [](std::vector<sparse_matrix*> spMat)
    {
        for (auto sp : spMat) {
            sp->setZero();
            sp->data().squeeze();
        }
    };

    // Symmetric gradient: 2*(0.5*(nabla u + nabla u^T), nabla v). The gradient
    // matrix is now
    //
    //     [A11 A12 A13]
    // K = [A21 A22 A23]
    //     [A31 A32 A33]
    //
    // where K21 = K12', K31 = K13', K32 = K23'.

    // Assemble diagonal matrices, and upper triangle off-diagonals that will
    // be recycled for the lower triangle of the full system matrix K.

    sparse_matrix A11(nd,nd);
    auto wA11 = WEAK_FORM( (2.0*ux.cwiseProduct(vx) + 1.0*uy.cwiseProduct(vy) + 1.0*uz.cwiseProduct(vz)) * Ve )
    fem.assemble_bilin( mesh, wA11, mu, A11, uDeg_, uDeg_ );
    sparse_matrix A22(nd,nd);
    auto wA22 = WEAK_FORM( (1.0*ux.cwiseProduct(vx) + 2.0*uy.cwiseProduct(vy) + 1.0*uz.cwiseProduct(vz)) * Ve )
    fem.assemble_bilin( mesh, wA22, mu, A22, uDeg_, uDeg_ );
    sparse_matrix A33(nd,nd);
    auto wA33 = WEAK_FORM( (1.0*ux.cwiseProduct(vx) + 1.0*uy.cwiseProduct(vy) + 2.0*uz.cwiseProduct(vz)) * Ve )
    fem.assemble_bilin( mesh, wA33, mu, A33, uDeg_, uDeg_ );

    sparse_matrix A12(nd,nd);
    auto wA12 = WEAK_FORM( 1.0*ux.cwiseProduct(vy) * Ve )
    fem.assemble_bilin( mesh, wA12, mu, A12, uDeg_, uDeg_ );
    sparse_matrix A13(nd,nd);
    auto wA13 = WEAK_FORM( 1.0*ux.cwiseProduct(vz) * Ve )
    fem.assemble_bilin( mesh, wA13, mu, A13, uDeg_, uDeg_ );
    sparse_matrix A23(nd,nd);
    auto wA23 = WEAK_FORM( 1.0*uy.cwiseProduct(vz) * Ve )
    fem.assemble_bilin( mesh, wA23, mu, A23, uDeg_, uDeg_ );

    // Full system matrix with Neumann for both velocity and pressure:
    A_ = sparse_matrix(3*nd+np, 3*nd+np);

    // A quarter row of the full matrix A_. Due to limitations of sparse block 
    // operations, need to fill Ar as a columns matrix, then transpose for A_.
    sparse_matrix Ar(3*nd+np, nd);   

    // First quarter row:
    Ar.middleRows(0*nd, nd) = A11.transpose();
    Ar.middleRows(1*nd, nd) = A12.transpose();
    Ar.middleRows(2*nd, nd) = A13.transpose();
    Ar.middleRows(3*nd, np) = -Bx_.transpose();
    A_.middleRows(0*nd, nd) = Ar.transpose();       
    freeSparseMemory_( {&A11} );       // free unused memory

    // Second row:
    Ar.middleRows(0*nd, nd) = A12;  // i.e., A21.transpose()
    Ar.middleRows(1*nd, nd) = A22.transpose();
    Ar.middleRows(2*nd, nd) = A23.transpose();
    Ar.middleRows(3*nd, np) = -By_.transpose();
    A_.middleRows(1*nd, nd) = Ar.transpose();       
    freeSparseMemory_( {&A12, &A22} );

    // Third row:
    Ar.middleRows(0*nd, nd) = A13;
    Ar.middleRows(1*nd, nd) = A23;
    Ar.middleRows(2*nd, nd) = A33.transpose();
    Ar.middleRows(3*nd, np) = -Bz_.transpose();
    A_.middleRows(2*nd, nd) = Ar.transpose();       
    freeSparseMemory_( {&A13, &A23, &A33} );  

    // Fourth row:
    Ar.middleRows(0*nd, nd) = Bx_;
    Ar.middleRows(1*nd, nd) = By_;
    Ar.middleRows(2*nd, nd) = Bz_;
    if (uDeg_ == 2) {      // P2-P1 element
        Ar.middleRows(3*nd, np) = par_["eps"]*M_.transpose();
    }
    else {              // P1-P1 stabilized
        sparse_matrix D(np, np);
        auto wD = WEAK_FORM( (h*h) * (ux.cwiseProduct(vx) + 
                                      uy.cwiseProduct(vy) + 
                                      uz.cwiseProduct(vz)) )
        fem.assemble_bilin( mesh, wD, {}, D );
        Ar.middleRows(3*nd, np) = (par_["alpha"]*D + par_["eps"]*M_).transpose();
        freeSparseMemory_( {&D} ); 
    }
    A_.middleRows(3*nd, np) = Ar.transpose();       
    freeSparseMemory_( {&Ar} ); 

    // Subset of K such that velocity has Dirichlet boundaries, pressure Neumann.
    subSparse(A_, sIdx_, sIdx_);    // uses tons of memory!

    std::cout << "  A: " << A_.rows() << " x " << A_.cols() << ", "
              << "number of nonzeros: " << A_.nonZeros() << std::endl;

    // Send the system matrix to the numerical solver. 
    // Stores reference of A_, not a copy!
    solver_.compute(A_);
}



int Flow::Stokes_( const LevelSet& lset, const VecXd& source, VecXd& flow )
{
    // To facilitate the convergence of the BiCGStab with variable viscosity, not
    // solving the Stokes if 1) variable viscosity, 2) RD iter > 0 (patterning 
    // requested), 3) growth factor concentration (source) zero.
    if (par_["Viscosity interior"] - par_["Viscosity exterior"] > FLT_MIN &&
        par_["RD iterations"] > 0.0 && source.cwiseAbs().maxCoeff() <= FLT_MIN)
        return 0;

    if (par_["Advection dt"] < FLT_MIN)
        return 0;

    chrono_time timeStart;
    tic( timeStart );

    //
    // Initialization, assemble matrices.
    // 

    // get viscosity distribution
    VecXd mu;
    Vec2f muMinMax;
    setViscosity( lset, par_, source, mu, muMinMax );
    std::cout << "  Viscosity min: " << muMinMax(0) << ", max: " << muMinMax(1) << std::endl;
    
    // fill the right-hand side
    VecXd FS = fillF_( lset, source, mu );
    std::cout << "  F: " << FS.rows() << " x " << FS.cols() << ". " 
              << "Mean: " << FS.mean() << ", max: " << FS.maxCoeff() << ", "
              << "min: " << FS.minCoeff() << ", hasNan(): " << FS.hasNaN() << std::endl;
    if (FS.norm() <= FLT_MIN) {
        std::cout << "  Load zero, nothing to solve." << std::endl;
        return 0;
    }
    // construct all the mutable matrices (left-hand side)
    constructMatrices_( lset, mu );
    std::cout << "  assembly took " <<  toc( timeStart ) << " ms." << std::endl;

  
    //
    // CPU solver, Eigen
    //
    int iter = 0;               // number of iterations the solver used
    double error = 0.0;         // residual error 
    tic( timeStart );

    std::cout << "  Solving (CPU)..." << std::flush;
    if (USE_INITIAL_GUESS)
        up_ = solver_.solveWithGuess(FS, up_);
    else
        up_ = solver_.solve(FS);

    iter = solver_.iterations();   
    error = solver_.error();

    std::cout << " finished after " << iter << " iterations and " 
              << toc(timeStart) << " ms."<< std::endl;
    std::cout << "  Error estimate: " << error << std::endl;
    size_t np = lset.getMesh().p().size();
    VecXd p = up_.tail(np);
    std::cout << "  Pressure mean " << p.mean() << ", max. " << p.maxCoeff() 
              << ", min. " << p.minCoeff() << std::endl;

    if (solver_.info() != Eigen::Success) {
        std::cerr << __FUNCTION__ << "(): Solving Stokes equation failed "
                  << "for unknown reason." << std::endl;
        return -1;
    }

    //
    // Copy interior node values (the non-zeros) to the full-domain vector.
    //
    size_t nd = lset.getMesh().getNDofs(uDeg_);
    VecXd up_full = VecXd::Zero(3*nd + np);
    for (size_t i=0; i<sIdx_.size(); i++) {
        up_full(sIdx_[i]) = up_(i);
    }
    for (size_t i=0; i<np; i++) {
        flow(i+0*np) = up_full(i);
        flow(i+1*np) = up_full(i+nd);
        flow(i+2*np) = up_full(i+2*nd);
        flow(i+3*np) = up_full(i+3*nd);
    }

    return 0;
}



void Flow::surfaceTension_( const float st, const int rho, const LevelSet& lset,
                            VecXd& fx, VecXd& fy, VecXd& fz )
{
    std::vector<float> k = lset.getInterfaceCurvature();

    #pragma omp parallel for
    for (size_t i=0; i<k.size(); i++) {
        auto& ife = lset.getInterfaceEdges();
        auto& phi = lset.getPhi();
        auto& mesh = lset.getMesh();

        auto e = mesh.edges( ife[i] );
        float t = phi[e(0)] / (phi[e(0)] - phi[e(1)]);

        auto& normals = lset.getNormals();
        auto& ifpArea = lset.getInterfaceNodeAreas();

        // The interface node at which k is given sits along the edge e; the values 
        // are distributed to the nodes of e.
        // Different conventions exist for the curvature sign. The curvatures in 
        // k are computed as \nabla \cdot n (as in Chen 2009), but need the reversed 
        // sign for the normal force (Brackbill 1992).
        Vec3f normal = normals[e(0)] / normals[e(0)].norm();
        Vec3f f = (1.0-t) * st * -k[i] * ifpArea[i] * rho * normal;
        fx(e(0)) += f(0);
        fy(e(0)) += f(1);
        fz(e(0)) += f(2);

        normal = normals[e(1)] / normals[e(1)].norm();
        f = t * st * -k[i] * ifpArea[i] * rho * normal;
        fx(e(1)) += f(0);
        fy(e(1)) += f(1);
        fz(e(1)) += f(2);
    }
}



VecXd Flow::fillF_( const LevelSet& lset, const VecXd& source, const VecXd& mu )
{
    auto& mesh = lset.getMesh();
    size_t np = mesh.p().size();        // number of P1 dofs
    size_t nd = mesh.getNDofs(uDeg_);   // number of dofs element (P1 or P2)

    // Load along each axis.
    VecXd fx = VecXd::Zero(np);
    VecXd fy = VecXd::Zero(np);
    VecXd fz = VecXd::Zero(np);

    //
    // Assign interfacial tension.
    //
    float st = par_["Surface tension"];
    if (st > 0.0f) {
        int rho = (int)par_["Node density"];
        surfaceTension_( st, rho, lset, fx, fy, fz );
    }

    //
    // Sources
    // 
    VecXd s = VecXd::Zero(np);
    for (int i : lset.getInteriorNodes()) {
        if (source(i) > par_["Growth threshold"])
            s(i) = (source(i) - par_["Level sink"]) * par_["GF growth"];
    }

    // Zero-mean sources, as required by the divergence theorem.
    // TODO: Should do actual integrations, now just assuming equal elements.
    s -= s.mean() * VecXd::Ones(np);

    // 
    // Fill the load-source vector.
    // TODO: Could save time by not reassembling matrices every iteration 
    // when using constant viscosity.
    // 
    auto& fem = lset.getAssembler();
    sparse_matrix Cx(nd,np), Cy(nd,np), Cz(nd,np);
    auto wCx = WEAK_FORM( Ve * u.cwiseProduct(vx) )
    fem.assemble_bilin( mesh, wCx, 2.0/3*mu, Cx, uDeg_, 1 );
    auto wCy = WEAK_FORM( Ve * u.cwiseProduct(vy) )
    fem.assemble_bilin( mesh, wCy, 2.0/3*mu, Cy, uDeg_, 1 );
    auto wCz = WEAK_FORM( Ve * u.cwiseProduct(vz) )
    fem.assemble_bilin( mesh, wCz, 2.0/3*mu, Cz, uDeg_, 1 );

    VecXd F = VecXd::Zero(3*nd+np);
    F.segment(0*nd, nd) = N_*fx + Cx*s;
    F.segment(1*nd, nd) = N_*fy + Cy*s;
    F.segment(2*nd, nd) = N_*fz + Cz*s;
    F.segment(3*nd, np) = M_*s;

    // Subvector of F according to the boundary conditions.
    VecXd Fd = subDense(F, sIdx_);

    return Fd;
}
