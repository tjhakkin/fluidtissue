//
// Viscosity profiles for growth; describes how the mass viscosity depends 
// on the growth factor concentration.
//

#pragma once

#include "fem/LevelSet.hpp"
#include "fem/MathTypes.hpp"


namespace FEM {

//
// Inverse profile is used in manuscript Figs. 4-5., direct profiles in Fig. 6.
// Just for testing, nothing fancy behind these. The profile is indicated by
// parameter 'Viscosity profile'.
//
enum ViscosityProfile { INVERSE, DIRECT, DIRECT_RATE4X };


//
// Fills element viscosity 'mu' based on morphogen concentration 'gf' and
// requested viscosity profile. 'muMinMax' stores the min. & max. mu values over
// the domain.
//
void setViscosity( const LevelSet& lset, const Parameters& par, const VecXd& gf, 
                   VecXd& mu, Vec2f& muMinMax )
{
    auto& mesh = lset.getMesh();

    double mu_h = par["Viscosity interior"];
    double mu_l = par["Viscosity exterior"];
    // exterior viscosity assumed to be constant:
    mu = mu_l * VecXd::Ones(mesh.t().size());

    double maxGF = 1.0;     // max. growth factor level
    if (par["RD iterations"] > 0.0) {
        // growth factor production equals the maximum concentration level
        maxGF = par["GF production"];
    }

    muMinMax = {FLT_MAX, FLT_MIN};
    for (int i : lset.getInteriorElements()) { 
        // element concentrations
        auto& t = mesh.t(i);
        Vec4d ec( gf(t(0)), gf(t(1)), gf(t(2)), gf(t(3)) );
        double c3 = ec.maxCoeff();

        // Element viscosity is determined by the maximum node concentration
        // in the element. Using element mean leads to the undesired effect 
        // of having very high viscosity at growth sources when Visc. int. large.
        if (par["Viscosity profile"] == INVERSE) {
            mu(i) = mu_l*(c3/maxGF) + mu_h*(1.0 - c3/maxGF);
        }
        if (par["Viscosity profile"] == DIRECT) {
            mu(i) = mu_h*(c3/maxGF) + mu_l*(1.0 - c3/maxGF);
        }
        if (par["Viscosity profile"] == DIRECT_RATE4X) {
            c3 *= 4.0;
            if (c3 > maxGF)
                c3 = maxGF;
            mu(i) = mu_h*(c3/maxGF) + mu_l*(1.0 - c3/maxGF);
        }

        if (muMinMax(0) > mu(i))
            muMinMax(0) = mu(i);
        if (muMinMax(1) < mu(i))
            muMinMax(1) = mu(i);
    }

}

}   // END namespace
