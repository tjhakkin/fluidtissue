//
// Initial level set shapes.
//
// Note: phi returned by each method may not have the correct distances over 
// the whole domain; only the zero is guaranteed to have the correct position.
//

#pragma once

#include "TetraMesh.hpp"
#include "MathTypes.hpp"


namespace FEM {

enum { ELLIPSOID, CUBOID, LEVEQUE };


//
// Sets initial level set as an origin centered ellipsoid.
//
void Ellipsoid( const Mesh3D& mesh, const float a, const float b, const float c, 
                VecXd& phi )
{
    for (size_t i=0; i<mesh.p().size(); i++) {
        auto& p = mesh.p(i);
        Vec3f pe = p.cwiseQuotient( Vec3f(a, b, c) );
        // Gives the correct sign but incorrect distance:
        phi(i) = pe.norm() - 1.0f;
    }
}

//
// Rectangular cuboid (currently cube with side length 'width' only!)
// TODO: an actual cuboid.
//
void Cuboid( const Mesh3D& mesh, const float width, const float height, const float depth, 
             VecXd& phi )
{
    for (size_t i=0; i<mesh.p().size(); i++) {
        auto& p = mesh.p(i).cwiseAbs();
        // temporary: use width only
        Vec3f d = p - Vec3f(width, width, width);
        if (d.maxCoeff() < 0.0f)
            phi(i) = p.maxCoeff() - width;
        else if (d.maxCoeff() > 0.0f)
            phi(i) = p.maxCoeff() - width;
        else
            phi(i) = 0.0;
    }
}

//
// Sphere of radius 0.15 centered at -(0.15, 0.15, 0.15) for the Enright's test.
// Note that the original Enright's test uses the unit box, hence the sphere
// was centered at (0.35, 0.35, 0.35).
//
void LeVeque( const Mesh3D& mesh, VecXd& phi )
{
    double rInv = 1 / 0.15;
    for (size_t i=0; i<mesh.p().size(); i++) {
        auto& p = (mesh.p(i) + Vec3f({0.15, 0.15, 0.15})) * rInv;
        // Gives the correct sign but incorrect distance:
        phi(i) = p.norm() - 1.0f;
    }
}

}   // END namespace
