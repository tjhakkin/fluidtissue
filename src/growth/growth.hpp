//
// Growth simulation (Stokes + RD) main program.
//

#pragma once

#include "Knots.hpp"
#include "Flow.hpp"
#include "RD.hpp"
#include "fem/Assembler.hpp"
#include "fem/LevelSet.hpp"


class Growth
{

public:
    Growth() : 
        m_id(0)    
    {}

    //
    // Initializes the simulations with given parameters and run ID.
    //
    int init( const Parameters& parameters, const int id );

    //
    // Executes the main simulation loop (Algorithm 1 in manuscript).
    //
    int run( const int iterations );


private:

    // Writes current simulation state at iteration i into output files.
    int storeModelState_( const int i );


    Knots               m_knots;        // differentiation, signalign centers
    LevelSet            m_lset;         // interface
    Flow                m_stokes;       // Stokes flow
    RD                  m_rd;           // patterning, growth factor

    FEM::VecXd          m_flow;         // Stokes velocity, pressure solution
    FEM::VecXd          m_activator;    // activator concentration in RD
    FEM::VecXd          m_inhibitor;    // inhibitor concentration in RD
    FEM::VecXd          m_growthFactor; // growth factor concentration in RD

    int                 m_id;           // run ID for file names
    Parameters          m_parameters;   // parameters read from the input file

};
