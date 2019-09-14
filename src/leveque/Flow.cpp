#include <iostream>

#include "Flow.hpp"
#include "common/Utils.hpp"

using namespace FEM;


void Flow::solve( const LevelSet& lset, const Parameters& par, const int iter,
                  VecXd& flow )
{
    LeVequeField_( lset, par, iter, flow );
}


//
//  PRIVATE METHODS
//


void Flow::LeVequeField_( const LevelSet& lset, const Parameters& par, 
                          const int iter, VecXd& flow )
{
    //
    // Generate the static field.
    //
    auto& mesh = lset.getMesh();
    size_t n = mesh.p().size();

    #pragma omp parallel for num_threads( nThreads_ )
    for (size_t i=0; i<n; i++) {
        // LeVeque's equations assume a unit cube centered at (0.5, 0.5, 0.5),
        // whereas our computational domain is centered at the origin, hence
        // need to translate:
        auto p = mesh.p(i).cast<double>() + Vec3d({0.5, 0.5, 0.5});
        
        flow(i)     = 2.0 * std::sin(M_PI*p(0)) * std::sin(M_PI*p(0)) * 
                      std::sin(2.0*M_PI*p(1)) * std::sin(2.0*M_PI*p(2));
        flow(i+n)   = -1.0 * std::sin(M_PI*p(1)) * std::sin(M_PI*p(1)) *
                      std::sin(2.0*M_PI*p(0)) * std::sin(2.0*M_PI*p(2));
        flow(i+2*n) = -1.0 * std::sin(M_PI*p(2)) * std::sin(M_PI*p(2)) *
                      std::sin(2.0*M_PI*p(0)) * std::sin(2.0*M_PI*p(1));
    }


    //
    // Time dependency.
    //
    int i = iter;
    double dt = par["Advection dt"];
    double T = 3.0;     // total time

    // Compute the velocity at the midpoint directly.
    double t = (i + 0.5) * dt;
    double g = std::cos(M_PI * t / T);
    flow *= g;
    
/*
    // Second order extrapolation for midpoint method: At iteration i returns
    // the velocity at time (i+1/2)*dt. Velocity is extrapolated as 
    //
    // u(i+1/2) = 3/2*u(i) - 1/2*u(i-1)
    //
    // with u(-1) = u(0). Could use true midpoint value here, but doing the 
    // explolation for consistency with other parts (surface tension, growth).
    double t0 = (i-1) * dt;
    if (t0 < 0.0)
        t0 = 0.0;
    double t1 = i * dt;
    double g0 = std::cos(M_PI * t0 / T);
    double g1 = std::cos(M_PI * t1 / T);
    flow = 0.5 * (3.0*flow*g1 - flow*g0);
*/
}
