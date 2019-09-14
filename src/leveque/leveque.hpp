//
// LeVeque/Enright's test main program.
//

#pragma once

#include "Flow.hpp"
#include "fem/Assembler.hpp"
#include "fem/LevelSet.hpp"

using namespace FEM;


class LeVeque
{

public:
    LeVeque() : 
        m_id(0)
    {}

    //
    // Initializes the simulations with given parameters and run ID.
    //
    int init( const Parameters& parameters, const int id );

    //
    // Executes the main simulation loop.
    //
    int run( const int iterations );


private:

    // Writes current simulation state at iteration i into output files.
    int storeModelState_( const int step );


    Flow                m_flow;         // LeVeque velocity field implementation
    LevelSet            m_lset;         // interface

    FEM::VecXd          m_flowField;    // local copy of the velocity field

    int                 m_id;           // run ID for file names
    Parameters          m_parameters;   // parameters read from the input file

};
