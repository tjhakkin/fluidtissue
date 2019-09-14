//
// The tetrahedral element. Implements the basis/shape functions and their 
// derivatives for polynomial degrees 1 and 2.
//

#pragma once

#include <functional>
#include <vector>
#include "MathTypes.hpp"

#define FEM_SHAPE( handle ) [](var_type x, var_type y, var_type z) { return handle; }


namespace FEM {

typedef float var_type;        // variable type for element functions
typedef std::function<var_type(var_type x, var_type y, var_type z)> fhandle;

struct Tetrahedron {
    // Initialize tethedral element with shape function polynomial degree 'deg'.
    Tetrahedron(int deg = 1) :
        degree_     (deg),
        N_          (0),
        edgeDofs_   (0)
    {
        if (degree_ == 1) {
            nDofs_ = 4;

            N_.resize(nDofs_);
            N_ = { FEM_SHAPE(1.0 - x - y - z),      // node 0 (0,0,0)
                   FEM_SHAPE(x),                    // node 1 (1,0,0)
                   FEM_SHAPE(y),                    // node 2 (0,1,0)
                   FEM_SHAPE(z) };                  // node 3 (0,0,1)

            dNdx_.resize(nDofs_);
            dNdx_ = { FEM_SHAPE(-1.0),
                      FEM_SHAPE(1.0),
                      FEM_SHAPE(0.0),
                      FEM_SHAPE(0.0) };
            dNdy_.resize(nDofs_);
            dNdy_ = { FEM_SHAPE(-1.0),
                      FEM_SHAPE(0.0),
                      FEM_SHAPE(1.0),
                      FEM_SHAPE(0.0) };
            dNdz_.resize(nDofs_);
            dNdz_ = { FEM_SHAPE(-1.0),
                      FEM_SHAPE(0.0),
                      FEM_SHAPE(0.0),
                      FEM_SHAPE(1.0) };
        }
        else if (degree_ == 2) {
            nDofs_ = 10;

            // edge bisector node pairs
            edgeDofs_ = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };

            N_.resize(nDofs_);
            N_ = { FEM_SHAPE( (1.0-x-y-z) * (1.0-2.0*x-2.0*y-2.0*z) ),  // node 0 (0,0,0)
                   FEM_SHAPE( x*(2.0*x-1.0) ),                  // node 1 (1,0,0)
                   FEM_SHAPE( y*(2.0*y-1.0) ),                  // node 2 (0,1,0)
                   FEM_SHAPE( z*(2.0*z-1.0) ),                  // node 3 (0,0,1)
                   FEM_SHAPE( 4.0*x*(1.0-x-y-z) ),              // edge (0,1)
                   FEM_SHAPE( 4.0*y*(1.0-x-y-z) ),              // edge (0,2)
                   FEM_SHAPE( 4.0*z*(1.0-x-y-z) ),              // edge (0,3)
                   FEM_SHAPE( 4.0*x*y ),                        // edge (1,2)
                   FEM_SHAPE( 4.0*x*z ),                        // edge (1,3)
                   FEM_SHAPE( 4.0*y*z ) };                      // edge (2,3)

            dNdx_.resize(nDofs_);
            dNdx_ = { FEM_SHAPE( 4.0*x + 4.0*y + 4.0*z - 3.0 ),
                      FEM_SHAPE( 4.0*x - 1.0 ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( -8.0*x - 4.0*y - 4.0*z + 4.0 ),
                      FEM_SHAPE( -4.0*y ),
                      FEM_SHAPE( -4.0*z ),
                      FEM_SHAPE( 4.0*y ),
                      FEM_SHAPE( 4.0*z ),
                      FEM_SHAPE( 0.0 ) };
            dNdy_.resize(nDofs_);
            dNdy_ = { FEM_SHAPE( 4.0*x + 4.0*y + 4.0*z - 3.0 ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( 4.0*y - 1.0 ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( -4.0*x ),
                      FEM_SHAPE( -4.0*x - 8.0*y - 4.0*z + 4.0 ),
                      FEM_SHAPE( -4.0*z ),
                      FEM_SHAPE( 4.0*x ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( 4.0*z ) };
            dNdz_.resize(nDofs_);
            dNdz_ = { FEM_SHAPE( 4.0*x + 4.0*y + 4.0*z - 3.0 ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( 4.0*z - 1.0 ),
                      FEM_SHAPE( -4.0*x ),
                      FEM_SHAPE( -4.0*y ),
                      FEM_SHAPE( -4.0*x - 4.0*y - 8.0*z + 4.0 ),
                      FEM_SHAPE( 0.0 ),
                      FEM_SHAPE( 4.0*x ),
                      FEM_SHAPE( 4.0*y ) };
            }
        else {
            return;
        }
    }


    // Returns handles to shape functions.
    const std::vector<fhandle>& N() const           { return N_; }

    // Returns handles to partial derivatives.
    const std::vector<fhandle>& dNdx() const        { return dNdx_; }
    const std::vector<fhandle>& dNdy() const        { return dNdy_; }
    const std::vector<fhandle>& dNdz() const        { return dNdz_; }

    const std::vector<Vec2i> getEdgeDofs() const    { return edgeDofs_; }


private:
    // shape function polynomial degree
    int                     degree_;
    
    // number of reference element degrees of freedom (dof)
    int                     nDofs_;

    // shape functions
    std::vector<fhandle>    N_;

    // partial derivatives
    std::vector<fhandle>    dNdx_;
    std::vector<fhandle>    dNdy_;
    std::vector<fhandle>    dNdz_;

    // edge dofs as edge node pairs (only for degree 2)
    std::vector<Vec2i>      edgeDofs_;

};

}   // END namespace
