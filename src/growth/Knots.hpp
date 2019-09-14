//
// Implements signaling centers that produce the growth factor.
//
// Terminology: 
// knot: a single position that has surpassed the threshold morphogen 
//       concentration for differentiation. 
// knot cluster: a collection of knots that have been clustered together.
// signaling center: a point object on the interface obtained from a knot cluster.
//
// Each signaling center gives rise to sources within the tissue interior near
// the signaling center.
//

#pragma once

#include "fem/MathTypes.hpp"
#include "fem/LevelSet.hpp"
#include "common/Parameters.hpp"


class Knots {

public:
    Knots() : 
        hK_(0.0f)
    {}

    //
    // Initialize variables, compute mesh parameter hK_.
    //
    void    init( const Parameters& par );

    //
    // Add new knots if morphogen level u surpasses a threshold (parameter 
    // 'Knots threshold') within hK_ distance from the interface. Individual knots 
    // are then clustered into signaling centers.
    //
    int     addKnots( const LevelSet& lset, const FEM::VecXd& u );

    //
    // Update signaling center positions with given velocity field.
    //
    int     updateSignalingCenters( const LevelSet& lset, const FEM::VecXd& flow );

    //
    // Generates vertices p and triangles t for visualizing the signaling centers as cubes.
    //
    void    knotsToMesh( std::vector<FEM::Vec3f>& p, std::vector<FEM::Vec3i>& t ) const;

    //
    // Returns the global indices of source nodes in the interior domain
    // induced by the signaling centers on the interface.
    //
    int     getSourceNodes( const LevelSet& lset, std::vector<int>& sourceNodes ) const;


private:

    // Computes the position of a mean knot representing each knot cluster.
    // Fills signalingCenters_.
    int    createSignalingCenters_( const LevelSet& lset,
                                    const std::vector<std::vector<int>>& clusters,
                                    const std::vector<FEM::Vec3f>& knots_p );

    // Fills 'cluster' with nodes belonging to the same cluster.
    // Removes the corresponding nodes from 'nodes'. Recursively calls itself.
    void    fillCluster_( const int cnode, const std::vector<FEM::Vec3f>& knots_p,
                          std::vector<int>& group, std::vector<int>& nodes, 
                          int depth ) const;

    // Combines spatially nearby knots into clusters. Calls fillCluster_().
    void    createKnotClusters_( const std::vector<FEM::Vec3f>& knots_p,
                                 std::vector<std::vector<int>>& clusters ) const;


    Parameters                      parameters_;        // model parameters

    float                           hK_;                // mesh parameter
    std::vector<FEM::Vec3f>         knots_p_;           // knot positions
    std::vector<std::vector<int>>   clusters_;          // knot clusters
    std::vector<FEM::Vec3f>         signalingCenters_;  // mean cluster knots

};
