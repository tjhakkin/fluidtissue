//
// Quadrature rules for the tetrahedral element.
//

#pragma once

#include <vector>
#include "MathTypes.hpp"
#include "Tetrahedron.hpp"


namespace FEM {

typedef double T;        // variable type for element functions

typedef Eigen::Matrix<T, Eigen::Dynamic, 1>    VecX;

struct Quadrature {

    Quadrature(int deg, const Tetrahedron& tet) :
        degree_     (deg),
        element_    (tet)
    {
        if (degree_ == 2) {
            nqp_ = 4;

            qp_ = {VecX(nqp_), VecX(nqp_), VecX(nqp_)};
            double a = (5.0 + 3.0*std::sqrt(5.0)) / 20.0;
            double b = (5.0 - std::sqrt(5.0)) / 20.0;
            qp_[0] << a, b, b, b;
            qp_[1] << b, a, b, b;
            qp_[2] << b, b, a, b;

            qw_ = VecX(nqp_);
            qw_ << 0.25, 0.25, 0.25, 0.25;
        }
        else if (degree_ == 3) {
            nqp_ = 10;

            qp_ = {VecX(nqp_), VecX(nqp_), VecX(nqp_)};
            qp_[0] << 0.7784952948213300, 0.0738349017262234, 0.0738349017262234, 0.0738349017262234, 0.4062443438840510, 0.4062443438840510, 0.4062443438840510, 0.0937556561159491, 0.0937556561159491, 0.0937556561159491;
            qp_[1] << 0.0738349017262234, 0.7784952948213300, 0.0738349017262234, 0.0738349017262234, 0.4062443438840510, 0.0937556561159491, 0.0937556561159491, 0.4062443438840510, 0.4062443438840510, 0.0937556561159491;
            qp_[2] << 0.0738349017262234, 0.0738349017262234, 0.7784952948213300, 0.0738349017262234, 0.0937556561159491, 0.4062443438840510, 0.0937556561159491, 0.4062443438840510, 0.0937556561159491, 0.4062443438840510;

            qw_ = VecX(nqp_);
            qw_ << 0.0476331348432089, 0.0476331348432089, 0.0476331348432089, 0.0476331348432089, 0.1349112434378610, 0.1349112434378610, 0.1349112434378610, 0.1349112434378610, 0.1349112434378610, 0.1349112434378610;
        }
        else if (degree_ == 4) {
            nqp_ = 20;

            qp_ = {VecX(nqp_), VecX(nqp_), VecX(nqp_)};
            qp_[0] << 0.9029422158182680, 0.0323525947272439, 0.0323525947272439, 0.0323525947272439, 0.2626825838877790, 0.6165965330619370, 0.2626825838877790, 0.6165965330619370, 0.2626825838877790, 0.6165965330619370, 0.0603604415251421, 0.0603604415251421, 0.0603604415251421, 0.0603604415251421, 0.0603604415251421, 0.0603604415251421, 0.3097693042728620, 0.3097693042728620, 0.3097693042728620, 0.0706920871814129;
            qp_[1] << 0.0323525947272439, 0.9029422158182680, 0.0323525947272439, 0.0323525947272439, 0.6165965330619370, 0.2626825838877790, 0.0603604415251421, 0.0603604415251421, 0.0603604415251421, 0.0603604415251421, 0.2626825838877790, 0.6165965330619370, 0.2626825838877790, 0.6165965330619370, 0.0603604415251421, 0.0603604415251421, 0.3097693042728620, 0.3097693042728620, 0.0706920871814129, 0.3097693042728620;
            qp_[2] << 0.0323525947272439, 0.0323525947272439, 0.9029422158182680, 0.0323525947272439, 0.0603604415251421, 0.0603604415251421, 0.6165965330619370, 0.2626825838877790, 0.0603604415251421, 0.0603604415251421, 0.6165965330619370, 0.2626825838877790, 0.0603604415251421, 0.0603604415251421, 0.2626825838877790, 0.6165965330619370, 0.3097693042728620, 0.0706920871814129, 0.3097693042728620, 0.3097693042728620;
        
            qw_ = VecX(nqp_);
            qw_ << 0.0070670747944695, 0.0070670747944695, 0.0070670747944695, 0.0070670747944695, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.0469986689718877, 0.1019369182898680, 0.1019369182898680, 0.1019369182898680, 0.1019369182898680;
        }
        else {}

        computeElementValues_();
    }


    // Returns quadrature points.
    const std::vector<VecX>& qp() const         { return qp_; }

    // Returns quadrature weights.
    const VecX& qw() const                      { return qw_; }

    // Returns shape function values at quadrature points.
    const std::vector<VecX>& phi() const        { return phi_; }

    // Returns shape function gradients at quadrature points.
    const std::vector<MatX>& gphi() const       { return gphi_; }


private:
    // shape function polynomial degree
    int                 degree_;

    Tetrahedron         element_;

    // number of quadrature points
    int                 nqp_;

    // quadrature points
    std::vector<VecX>   qp_;

    // quadrature weights
    VecX                qw_;

    // shape function values at quadrature points
    std::vector<VecX>   phi_;

    // shape function partial derivatives at quadrature points
    std::vector<MatX>   gphi_;


    // Compute values of shape functions and shape function partial derivatives
    // at the quadrature points.
    void computeElementValues_()
    {
        // Shape function values
        auto& N = element_.N();
        phi_.resize(N.size());

        for (size_t i=0; i<N.size(); i++) {
            phi_[i] = VecX(nqp_);
            for (int j=0; j<nqp_; j++)
                phi_[i](j) = N[i](qp_[0](j), qp_[1](j), qp_[2](j));
        }

        // Shape function partial derivative values.
        auto& dNdx = element_.dNdx();
        auto& dNdy = element_.dNdy();
        auto& dNdz = element_.dNdz();
        gphi_.resize(dNdx.size());

        for (size_t i=0; i<gphi_.size(); i++) {
            gphi_[i] = MatX(3, nqp_);           // dN/dx, dN/dy, dN/dz -> 3
            for (int j=0; j<nqp_; j++) {
                gphi_[i](0,j) = dNdx[i](qp_[0](j), qp_[1](j), qp_[2](j));
                gphi_[i](1,j) = dNdy[i](qp_[0](j), qp_[1](j), qp_[2](j));
                gphi_[i](2,j) = dNdz[i](qp_[0](j), qp_[1](j), qp_[2](j));
            }
        }

    }
};

}   // END namespace
