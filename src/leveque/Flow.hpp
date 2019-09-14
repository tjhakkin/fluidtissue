#pragma once

#include "fem/MathTypes.hpp"
#include "fem/Assembler.hpp"
#include "fem/LevelSet.hpp"
#include "common/Parameters.hpp"


class Flow {
 
public:
    Flow() : 
        nThreads_(std::thread::hardware_concurrency())
    {}

    //
    // Fills 'flow' with the velocity field. Calls LevequeField_(), nothing else.
    //
    void    solve( const LevelSet& lset, const Parameters& par, const int iter,
                   FEM::VecXd& flow );


private:

    // Fills 'flow' at time (iter+0.5)*dt. Called from solve().
    void    LeVequeField_( const LevelSet& lset, const Parameters& par, 
                           const int iter, FEM::VecXd& flow );

    int nThreads_;          // max. number of CPU cores to use

};
