//
// Solves the three-morphogen reaction-diffusion system to obtain patterning
// and the growth factor distribution.
//

#pragma once

#include "Knots.hpp"
#include "fem/LevelSet.hpp"
#include "fem/MathTypes.hpp"
#include "fem/TetraMesh.hpp"
#include "common/Parameters.hpp"


class RD {

public:
    RD() : 
        m_nodes(0)
    {}

    const FEM::VecXd&       getActivator() const    { return m_activator; }
    const FEM::VecXd&       getInhibitor() const    { return m_inhibitor; }
    const FEM::VecXd&       getGrowthFactor() const { return m_growthFactor; }
    const std::vector<int>& getNodes() const        { return m_nodes; }
    const FEM::Mesh3D&      getMesh() const         { return m_meshRD; } 

    //
    // Set initial morphogen concentrations.
    //
    void init( const LevelSet& lset, const Parameters& par );
    
    //
    // Solves reaction-diffusion for patterning and growth factor for growth.
    //
    int solve( const LevelSet& lset, const Knots& knots );


private:

    // Fills the Laplacian L = invM * K, where invM the mass inverse and K the
    // stiffness matrix. Fills global-to-local node index map 'imap'.
    int constructMatrices_( const LevelSet& lset, std::map<int,int>& imap,
                            FEM::sparse_matrix& L );


    Parameters          m_parameters;   // model parameters

    std::vector<int>    m_nodes;        // interior node indices (full mesh)
    FEM::Mesh3D         m_meshRD;       // interior elements mesh

    FEM::VecXd          m_activator;    // activator concentration in RD
    FEM::VecXd          m_inhibitor;    // inhibitor concentration in RD
    FEM::VecXd          m_growthFactor; // growth factor concentration in RD

};
