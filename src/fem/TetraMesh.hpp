//
// Representation of the computational mesh structure and the associated FEM
// dofs, basis/shape functions and their derivatives. 
//
// ** Basic usage **
//
// Create an empty mesh and then populate it
//
// Mesh3D mesh;
// mesh.createRegularMesh(pars, nThreads);
//
// where 'pars' is the parameters struct and 'nThreads' the maximum number of 
// threads used by the triangulator. See LevelSet::init() for an example. To 
// access basic mesh attributes call, e.g.,
//
// auto& vertices = mesh.p(); 
// auto& vertex = mesh.p(n);
//
// to retrieve a vector of all mesh vertices and the nth vertex, respectively. 
//

#pragma once

#include <cfloat>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>

#include "MathTypes.hpp"
#include "Tetrahedron.hpp"
#include "common/Utils.hpp"
#include "common/Parameters.hpp"


namespace FEM {

class Mesh3D {

public:

    Mesh3D() : 
        element_    {1,2}        // P1, P2 elements
    {}

    //
    // Retrieval methods for basic mesh attributes (vertices, tetrahedra, ...)
    //
    const auto&     p()             const   { return p_; }          // vertices
    const auto&     p(int i)        const   { return p_[i]; }
    const auto&     t()             const   { return t_; }          // tetrahedra
    const auto&     t(int i)        const   { return t_[i]; }       
    const auto&     edges()         const   { return edges_; }      // edges
    const auto&     edges(int i)    const   { return edges_[i]; }
    const auto&     p2t()           const   { return p2t_; }        // nodes to tetrahedra
    const auto&     p2t(int i)      const   { return p2t_[i]; }
    const auto&     t2e()           const   { return t2e_; }        // node to edges
    const auto&     t2e(int i)      const   { return t2e_[i]; }
    const auto&     b_tris()        const   { return b_tris_; }     // exterior boundary triangles
    const auto&     b_tris(int i)   const   { return b_tris_[i]; }
    const auto&     b_edgeIDs()     const   { return b_edgeIDs_; }  // exterior boundary edge IDs
    const auto&     b_edgeIDs(int i) const  { return b_edgeIDs_[i]; }

    //
    // Store methods for vertices, tetrahedra and node-to-tetradra.
    //
    void    set_p(std::vector<Vec3f>& p)                                    { p_ = p; }
    void    set_t(std::vector<Vec4i, Eigen::aligned_allocator<Vec4i>>& t)   { t_ = t; }
    void    set_p2t(std::vector<VecXi>& p2t)                                { p2t_ = p2t; }

    //
    // Creates an origin-centered cube/box mesh of tetrahedral elements with 
    // regularly spaced nodes.
    // Box dimensions are determined by domain width, height, depth -parameters.
    // The total number of nodes is controlled by the node density parameter and
    // the requested dimensions.
    //
    int     createRegularMesh( const Parameters& parameters, int nThreads );

    //
    // Returns handles to shape functions of degree d.
    //
    const std::vector<fhandle>&     N(int d) const      { return element_[d-1].N(); }

    //
    // Returns handles to partial derivatives.
    //
    const std::vector<fhandle>&     dNdx(int d) const   { return element_[d-1].dNdx(); }
    const std::vector<fhandle>&     dNdy(int d) const   { return element_[d-1].dNdy(); }
    const std::vector<fhandle>&     dNdz(int d) const   { return element_[d-1].dNdz(); }

    //
    // Returns the reference element of given degree.
    //
    const Tetrahedron&      getElement(int d) const      { return element_[d-1]; }

    //
    // Returns the interior dof indices.
    //
    const std::vector<int>& getInteriorDofIndices(int d) const  { return interiorDofs_[d-1]; };

    //
    // Returns the number of global dofs for given element polynomial degree d.
    //
    int getNDofs(int d) const
    {
        if (d == 1)
            return p_.size();
        else if (d == 2)
            return p_.size() + edges_.size();
        else
            return -1;
    }

    //
    // Returns dof indices of element k.
    //
    inline const VecXi getDofIndices(int k) const   { return dofIndices_[k]; }

    //
    // Returns the approximat memory footprint of the mesh.
    //
    int sizeInMemory() const
    {
        size_t m = p_.size() * sizeof(decltype(p_)::value_type) + 
                   t_.size() * sizeof(decltype(t_)::value_type) + 
                   edges_.size() * sizeof(decltype(edges_)::value_type) + 
                   4 * t_.size() * sizeof(decltype(p2t_)::value_type) + 
                   4 * t_.size() * sizeof(decltype(t2e_)::value_type) + 
                   b_tris_.size() * sizeof(decltype(b_tris_)::value_type);
        return m / (1024*1024);
    }

    //
    // Given a position, recursively find the nodes within given radius from the 
    // position.
    // - if radius <= 0 returns the closest node only
    // - closest node always at index 0
    //
    inline std::vector<int> getNodesWithinRadius( const Vec3f& po, 
                                                  const float radius = 0.0f ) const
    {
        std::vector<int> nodes;
        float dmin = FLT_MAX;
        int imin = -1;

        for (size_t i=0; i<p_.size(); ++i) {
            float d = (p_[i] - po).norm();
            if (d < dmin) {
                dmin = d;
                imin = i;
            }
            if (d <= radius)
                nodes.push_back(i);
        }

        nodes.erase(std::remove(nodes.begin(), nodes.end(), imin), nodes.end());
        nodes.insert(nodes.begin(), imin);

        return nodes;
    }

    //
    // Given a point finds the closet mesh node and its connected neigbours.
    //
    inline void getClosestNode( const Vec3f& po, std::vector<int>& nodes ) const
    {
        float dmin = FLT_MAX;
        int imin = -1;

        for (size_t i=0; i<p_.size(); ++i) {
            float d = (p_[i] - po).norm();
            if (d < dmin) {
                dmin = d;
                imin = i;
            }
        }
        nodes.push_back( imin );
        
        // Add the neighbors connected to imin.
        // TODO: On a large mesh this might be inefficient.
        for (auto& e : edges_) {
            if (e(0) == imin)
                nodes.push_back( e(1) );
            if (e(1) == imin)
                nodes.push_back( e(0) );
        }
    }

    //
    // Retrieves tetrahedron barycentric coordinates of given point.
    //
    inline void getTetrahedronBC( const std::vector<Vec3f>& v, const Vec3f& po, 
                                  Vec4f& coeff ) const
    {
        Mat3f M;
        M << v[0]-v[3], v[1]-v[3], v[2]-v[3];

        Vec3f c = M.inverse() * (po - v[3]);
        coeff = { c(0), c(1), c(2), 1.0f-c.sum() };
    }

    //
    // Given point 'po', fills the barycentric coordinates 'bc' and tetrahedron 'tt'
    // for the enclosing tetrahedon.
    // Returns 0 if enclosing tetrahedron found, else -1.
    //
    inline int getEnclosingTetrahedron( const Vec3f& po, Vec4f& bc, Vec4i& tt ) const
    {
        tt = { -1, 0, 0, 0 };                // tetrahedron containing the knot
        bc = { 0.0f, 0.0f, 0.0f, 0.0f };    // barycentric coefficients

        // Get the mesh node closest to the knot and its neighbours.
        std::vector<int> nodes;
        getClosestNode( po, nodes );

        for (int node : nodes ) {
            VecXi tetras = p2t_[ node ];

            for (int j=0; j<tetras.size(); ++j) {
                tt = t_[ tetras(j) ];

                std::vector<Vec3f> tp = { p_[tt(0)], p_[tt(1)], p_[tt(2)], p_[tt(3)] };
                getTetrahedronBC( tp, po, bc );

                if ( bc.maxCoeff() <= 1.0f && bc.minCoeff() >= 0.0f )
                    return 0;
            }
        }

        return -1;
    }


private:

    // Performs Delaunay triangulation on mesh vertices, fill the mesh fields.
    int Delaunay_3D( const Parameters& parameters, const int nThreads );

    // Given node density and domain width, height, depth parameters, computes the
    // number of mesh nodes along each axis and the final origin-centered mesh
    // dimensions (width, height, depth).
    void initMesh( const Parameters& parameters, Vec3i& nodes, Vec3f& dimensions );

    // Initiazes mesh dofs (degrees of freedom).
    void initDofs();


    // mesh attributes
    std::vector<Vec3f>                                  p_;         // vertices
    std::vector<Vec4i, Eigen::aligned_allocator<Vec4i>> t_;         // tetrahedra
    std::vector<Vec2i>                                  edges_;     // tetrahedra edges
    std::vector<VecXi>                                  p2t_;       // node to tetrahedra
    std::vector<Vec6i>                                  t2e_;       // tetrahedra to edges
    std::vector<Vec3i>                                  b_tris_;    // exterior boundary triangles
    std::vector<Vec3i>                                  b_edgeIDs_; // exterior boundary edge IDs

    Tetrahedron             element_[2];        // reference elements
    std::vector<VecXi>      dofIndices_;        // tetrahedron -> dofs maps
    std::vector<int>        interiorDofs_[2];   // list of interior/non-boundary dof indices

};

}   // END namespace
