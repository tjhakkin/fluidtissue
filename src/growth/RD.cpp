#include "RD.hpp"
#include "fem/Assembler.hpp"
#include "common/Utils.hpp"

using namespace FEM;


void RD::init( const LevelSet& lset, const Parameters& par )
{
    m_parameters = par;

    auto& mesh = lset.getMesh();
    size_t np = mesh.p().size();

    m_activator = VecXd::Zero(np);
    m_inhibitor = VecXd::Zero(np);
    m_growthFactor = VecXd::Zero(np);

    // Initial activator as a linear ramp.
    // TODO: Not supposed to be hard-coded here!
    for (auto i : lset.getInteriorNodes())
        m_activator(i) = mesh.p(i)(0) + 1.0;
}



int RD::solve( const LevelSet& lset, const Knots& knots )
{
    int nIter = (int)m_parameters["RD iterations"];
    double Dt = m_parameters["RD dt"];
    if (nIter == 0 || Dt < FLT_MIN)
        return 0;

    sparse_matrix  L;                   // Laplacian as invM * K
    std::map<int,int> iMap;             // global-to-local node map
    if (constructMatrices_(lset, iMap, L))
        return -1;

    // Copy concentrations to a local vector restricted to the inner domain.
    size_t np = L.rows();               // number of inner nodes, also L.cols().
    VecXd u(np), v(np), g(np);
    for (auto const& k : iMap) {
        u(k.second) = m_activator(k.first);
        v(k.second) = m_inhibitor(k.first);
        g(k.second) = m_growthFactor(k.first);
    }

    // To save time, consider activator-inhibitor patterning inactive when the
    // concentrations near zero.
    // NOTE: Should be disabled if background sources (actBG) enabled!
    bool RD_alive = true;
    if (m_activator.cwiseAbs().maxCoeff() < FLT_MIN && m_inhibitor.cwiseAbs().maxCoeff() < FLT_MIN)
        RD_alive = false;

    //
    // Growth sources are implemented by first setting the corresponding mesh
    // nodes to value 'GF production', then solving growth factor diffusion over 
    // the whole inner domain as usual, then setting the source values again.
    //

    // Get source node indices.
    std::vector<int> srcNodes; 
    knots.getSourceNodes( lset, srcNodes );
    for (size_t i=0; i<srcNodes.size(); i++)
        srcNodes[i] = iMap[srcNodes[i]];        // convert to local indexing

    // Initialize the source node values.
    for (int j : srcNodes)
        g(j) = m_parameters["GF production"];

    //
    // Solve the equations with explicit Euler method for dt*nIter time.
    //
    for (int i=0; i<nIter; i++) {
        auto& p = m_parameters;

        // Activator, inhibitor
        if (RD_alive) {
            VecXd u2 = u.cwiseProduct(u);                                           // u^2                
            VecXd de = VecXd::Ones(v.rows()) + p["Inh"]*v + p["GF inhibition"]*g;   // 1 + beta12*v + beta13*g
            VecXd un = u + Dt*( -p["Da"]*L*u + p["Act"]*u2.cwiseQuotient(de) - p["DegA"]*u );
            VecXd vn = v + Dt*( -p["Di"]*L*v + u2 - p["DegI"]*v );
            u = un;
            v = vn;
        }

        // Growth factor
        g = g + Dt*( -p["GF diffusion"]*L*g - p["GF degradation"]*g );
        for (int j : srcNodes)
            g(j) = p["GF production"];
    }

    // Copy local concentrations back to the full domain vector.
    m_activator.setZero();
    m_inhibitor.setZero();
    m_growthFactor.setZero();
    for (auto const& k : iMap) {
        m_activator(k.first) = u(k.second);
        m_inhibitor(k.first) = v(k.second);
        m_growthFactor(k.first) = g(k.second);
    }

    std::cout << "  Max. act: " << m_activator.maxCoeff() 
              << ", inh: " << m_inhibitor.maxCoeff()
              << ", gf: " << m_growthFactor.maxCoeff() << std::endl;
    
    return 0;
}


//
//  PRIVATE METHODS
//


int RD::constructMatrices_( const LevelSet& lset, std::map<int,int>& imap, 
                            sparse_matrix& L )
{
    //
    // Construct the interior domain submesh.
    //

    auto& mesh = lset.getMesh();
    m_nodes = lset.getInteriorNodes();

    // Node map from global node index to elements-local index.
    imap.clear();
    for (size_t i=0; i<m_nodes.size(); i++)
        imap[ m_nodes[i] ] = i;

    // Set elements with the new node indexing.
    std::vector<int> elements = lset.getInteriorElements();
    std::vector<Vec4i, Eigen::aligned_allocator<Vec4i>> tn( elements.size() );
    for (size_t i=0; i<elements.size(); i++) {
        auto& t = mesh.t( elements[i] );
        tn[i] = { imap[t(0)], imap[t(1)], imap[t(2)], imap[t(3)] };
    }

    // Update RD mesh elements and nodes.
    // NOTE: m_meshRD.p2t not filled, hence no calling assembler grad_nodes().
    m_meshRD.set_t(tn);
    std::vector<Vec3f> pn( m_nodes.size() );
    for (size_t i=0; i<pn.size(); i++)
        pn[i] = mesh.p( m_nodes[i] );
    m_meshRD.set_p(pn);

    //
    // Assemble FEM matrices within the interior mesh, with Neumann boundaries.
    //

    Assembler assembler;
    size_t n = m_meshRD.p().size();
    // mass matrix, bilinear form (u, v)
    sparse_matrix M(n,n);
    auto wM = WEAK_FORM( u.cwiseProduct(v) )
    if (assembler.assemble_bilin( m_meshRD, wM, {}, M ))
        return -1;

    // stiffness matrix, bilinear form (nabla u, nabla v)
    sparse_matrix K(n,n);
    auto wK = WEAK_FORM( ux.cwiseProduct(vx) + uy.cwiseProduct(vy) + uz.cwiseProduct(vz) )
    if (assembler.assemble_bilin( m_meshRD, wK, {}, K ))
        return -1;

    // Approximate mass inverse by lumping mass:
    sparse_matrix invM;
    assembler.lump_mass( M, invM );

    L = invM * K;       // Laplacian.

    return 0;
}
