// 
// Solves the Stokes flow with variable viscosity, masss sources and body forces
// using either the stable P2-P1 element, or the stabilized Brezzi-Pitkäranta
// P1-P1 element.
//

#pragma once

#include "fem/MathTypes.hpp"
#include "fem/Assembler.hpp"
#include "fem/LevelSet.hpp"
#include "common/Parameters.hpp"


class Flow {
 
public:
    Flow() : 
        uDeg_(1)
    {}

    //
    // Assembles the immutable FEM matrices for solving the Stokes flow,
    // initializes the numerical solver.
    //
    int     init( const LevelSet& lset, const Parameters& par );

    //
    // Solves the Stokes flow for velocity and pressure with given source field 
    // indicating growth sources.
    //
    int     solve( const LevelSet& lset, const FEM::VecXd& source, FEM::VecXd& flow );


private:

    // Assembles the mutable left-hand side Stokes system matrices.
    void    constructMatrices_( const LevelSet& lset, const FEM::VecXd& source );
 
    // The solver, called by solve().
    int     Stokes_( const LevelSet& lset, const FEM::VecXd& source, FEM::VecXd& flow );

    // Computes surface tension at the interface nodes, fills force components f*.
    void    surfaceTension_( const float st, const int rho, const LevelSet& lset,
                             FEM::VecXd& fx, FEM::VecXd& fy, FEM::VecXd& fz );

    // Returns the right-hand side of the Stokes system.
    FEM::VecXd  fillF_( const LevelSet& lset, const FEM::VecXd& source, const FEM::VecXd& mu );


    Parameters              par_;           // model parameters

    std::vector<int>        sIdx_;          // indices of solution nodes
    int                     uDeg_;          // velocity space polynomial degree

    FEM::sparse_matrix      A_;             // Stokes system matrix; 4n x 4n
    FEM::sparse_matrix      M_, N_;         // (u,v); n x n (P1 stab., P2)
    FEM::sparse_matrix      Bx_, By_, Bz_;  // (u,vx), (u,vy), (u,vz); n x n
    FEM::VecXd              up_;            // velocity-pressure -solution

    // Solver, no preconditioner:
    Eigen::BiCGSTAB<FEM::sparse_matrix, Eigen::IdentityPreconditioner> solver_;
};
