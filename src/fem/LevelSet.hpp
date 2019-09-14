//
// A simple level set implementation on a tetrahedral mesh.
//
// ** Basic usage **
//
// Given parameters object 'par' indicating mesh parameters, initial level set
// configuration, initialize with
// 
// LevelSet lset;
// lset.init(par);
//
// The mesh and the current level set values at the mesh nodes can be obtained by
//
// Mesh3D mesh = lset.getMesh();
// VecXd phi = lset.getPhi();
//
// Given velocity field 'u', the interface position is updated by calling
//
// lset.update(u);
//
// which, after advecting the interface, recomputes the interface/mesh distances 
// to keep 'phi' as a distance function. To reduce mass loss caused by distance
// recomputation, the final values of 'phi' at the interface-crossing edge nodes 
// is mixed with non-reset values using weights indicated by PHI_RESET_WEIGHTS.
//

#pragma once

#include "Assembler.hpp"
#include "EuclideanDistance.hpp"
#include "TetraMesh.hpp"
#include "common/Parameters.hpp"


class LevelSet {

public:
    LevelSet() : 
        nThreads_(std::thread::hardware_concurrency())     
    {}

    //
    // Initializes the computational domain, assembles immutable FEM matrices.
    //
    int     init( const Parameters& pars );

    //
    // Solves the advection equation on phi_ to update the interface position.
    //
    int     update( const FEM::VecXd& flow );

    //
    // Returns interpolated level set (phi) and level set normal (normal) at point p.
    //
    int     interpolatePhi( const FEM::Vec3f& p, float& phi, FEM::Vec3f& normal ) const;

    //
    // Returns indices of elements in the negative level set.
    //
    std::vector<int>            getInteriorElements() const;

    //
    // Returns indices of nodes with negative level set.
    //
    std::vector<int>            getInteriorNodes() const;

    //
    // Returns the curvature of phi interpolated at the interface nodes as
    // k = \nabla \cdot (\nabla \phi / |\nabla \phi|)
    //
    std::vector<float>          getInterfaceCurvature() const;

    //
    // Returns the reconstructed interface surface mesh.
    //
    const FEM::Mesh3D&          getInterface() const            { return ifMesh_; }

    //
    // Returns interface node -> edge map, i.e., mesh edges that are crossed
    // by the interface.
    // TODO: This is a bit badly/misleadingly named.
    //
    const std::vector<int>&     getInterfaceEdges() const       { return ifp2e_; }

    //
    // Returns aggregate areas of triangles associated with each interface node.
    //
    const std::vector<float>&   getInterfaceNodeAreas() const   { return ifpArea_; }     

    //
    // Returns the current level set.
    //
    const FEM::VecXd&           getPhi() const                  { return phi_; }

    //
    // Returns level set normals.
    //
    const std::vector<FEM::Vec3f>&  getNormals() const          { return normals_; }

    //
    // Returns the mesh on which phi is computed.
    //
    const FEM::Mesh3D&          getMesh() const                 { return mesh_; }

    //
    // Returns FEM assembler.
    //
    const FEM::Assembler&       getAssembler() const            { return assembler_; }


private:

    // Find interface-mesh intersection points.
    void    findInterfaceNodes_();

    // Fills ifMesh_, ifElements_, ifp_area_, called after level set reset.
    int     reconstructInterface_();

    // Resets level set back to signed distance function under the assumption
    // that the correct interface position can be computed from the current
    // level set (signs must be correct).
    void    resetPhi_();

    Parameters              par_;           // model parameters

    FEM::Mesh3D             mesh_;          // mesh created at init()
    FEM::Assembler          assembler_;     // FEM matrix assembler
    FEM::VecXd              phi_;           // distances to interface, weighted
    FEM::VecXd              phi_ifn_;       // distances to interface, full reset
    FEM::sparse_matrix      M_;             // FEM mass matrix for solveAdvection()
    std::vector<FEM::Vec3f> normals_;       // level set normals

    std::vector<FEM::Vec3f> ifp_;           // interface node positions along edges
    std::map<int,int>       e2ifp_;         // edge ID -> interface node map
    std::vector<int>        ifp2e_;         // interface node -> edge map

    FEM::Mesh3D             ifMesh_;        // reconstructed interface surface
    std::vector<int>        ifElements_;    // interface crossing elements
    std::vector<float>      ifpArea_;       // total area of triangles associated with interface nodes

    int                     nThreads_;      // max. number CPU threads
    EuclideanDistance       dist_;          // computes the interface/mesh distances to reset phi_, phi_ifn

    // With Jacobi/diagonal preconditioner (default).
    Eigen::BiCGSTAB<FEM::sparse_matrix> solver_;
};
