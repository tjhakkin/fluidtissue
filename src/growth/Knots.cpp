#include <cfloat>

#include "Knots.hpp"
#include "fem/TetraMesh.hpp"

using namespace FEM;


namespace {
    // Knots that are within CLUSTER_DIAMETER * hK_ distance from eachother will form a cluster.
    // TODO: Need a less arbitrary criterion, something derived from the patterning dynamics.
    const double CLUSTER_DIAMETER = 4.0;

    // Average distance of a source node from the interface as a multiple of hK_.
    // To avoid numerical issues, value slightly larger than 1.0 probably a good choice.
    const double SOURCE_DISTANCE = 1.5;

    // Side length of the cube to visualize signaling centers as a multiple of hK_.
    const double CUBE_SIZE = 1.1;
}



// TODO: Currently mesh parameter hK_ equals the length of a diagonal mesh edge,
// i.e., assuming regular mesh. Should be e.g. the diameter of mean excircle.
void Knots::init( const Parameters& par )
{
    parameters_ = par;

    int density = (int)par["Node density"];             // nodes per unit cube
    float nMax = std::floor( std::cbrt(density) );      // ...per side
    if (nMax < FLT_MIN)
        nMax = 1.0;
    hK_ = 1.0f / nMax;
}



// Knot formation is currently limited by min./max. angle relative to the apical 
// direction. Ideally, the formation should be controlled by the patterning
// dynamic, rather than artificial parameters like these.
int Knots::addKnots( const LevelSet& lset, const VecXd& u )
{
    //
    // Add knots if concentration 'u' above a threshold with hK_ distance from
    // the interface & passes the angle test.
    //
    std::vector<Vec3f> knots_p;       // new knot positions
    for (int i=0; i<u.rows(); ++i) 
    {
        auto& par = parameters_;
        auto& phi = lset.getPhi();
        if (u(i) > par["Knots threshold"] && std::abs(phi(i)) < hK_ ) 
        {
            auto& mesh = lset.getMesh();
            auto& normals = lset.getNormals();
            Vec3f n = normals[i] / normals[i].norm();
            Vec3f p = mesh.p(i) - phi(i)*n;

            // By present convention, the apical direction is negative:
            Vec3f apicalDir = Vec3f(0.0f, 0.0f, -1.0f);
            float dp = std::acos( n.dot(apicalDir) );
            if (dp > par["Knots angle max"] || dp < par["Knots angle min"])
                continue;

            knots_p.push_back(p);
            std::cout << "  Added knot at " << p.transpose() << std::endl;
        }
    }

    if (knots_p.size() == 0)
        return 0;

    //
    // Construct knot clusters and the corresponding signaling centers.
    //
    std::vector<std::vector<int>> clusters;
    createKnotClusters_( knots_p, clusters );
    if (createSignalingCenters_( lset, clusters, knots_p ))
        return -1;

    knots_p_.insert( knots_p_.end(), knots_p.begin(), knots_p.end() );
    clusters_.insert( clusters_.end(), clusters.begin(), clusters.end() );

    std::cout << "  Number of knot clusters " << clusters_.size() << ", with (";
    for (size_t i=0; i<clusters_.size(); i++) {
        std::cout << clusters_[i].size();
        if (i < (clusters_.size()-1))
            std::cout << ", ";
    }
    std::cout << ") knots." << std::endl;

    return 0;
}



int Knots::updateSignalingCenters( const LevelSet& lset, const VecXd& flow )
{
    float dt = (float)parameters_["Advection dt"];
    if (dt < FLT_MIN)
        return 0;

    for (Vec3f& p : signalingCenters_) {
        auto& mesh = lset.getMesh();

        float lset_p;       // level set at p
        Vec3f normal_p;     // level set normal at p
        if (lset.interpolatePhi(p, lset_p, normal_p))
            return -1;
        Vec4f bc;           // barycentric coordinates of p inside t
        Vec4i t;            // p enclosing tetrahedron
        if (mesh.getEnclosingTetrahedron(p, bc, t))
            return -1;

        // Interpolate velocity at the signaling center.
        size_t n = mesh.p().size();
        Vec3f up( 0.0f, 0.0f, 0.0f );
        for (int j=0; j<4; j++)
            up(0) += bc(j) * flow(t(j));        // x component
        for (int j=0; j<4; j++)
            up(1) += bc(j) * flow(t(j)+n);      // y component
        for (int j=0; j<4; j++)
            up(2) += bc(j) * flow(t(j)+2*n);    // z component


        // Update signaling center position, normal component only:
        up = up.dot( normal_p ) * normal_p;
        p += dt * up;

        // Signaling centers may drift off the interface over time, correct 
        // the position if needed.
        if (lset.interpolatePhi(p, lset_p, normal_p))
            return -1;
        if (std::fabs(lset_p) > hK_) {
            p -= lset_p * normal_p;
            float lset_old = lset_p;
            if (lset.interpolatePhi(p, lset_p, normal_p))
                return -1;
            std::cout << "  Corrected knot position (old phi " << lset_old
                      << ", new phi " << lset_p << ")" << std::endl;
        }
    }

    return 0;
}



void Knots::knotsToMesh( std::vector<Vec3f>& p, std::vector<Vec3i>& t ) const
{
    Eigen::MatrixXf cubeVertices( 3, 8 );       // 8 vertices in a cube
    cubeVertices << -1.0f,  1.0f, 1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f,
                    -1.0f, -1.0f, 1.0f,  1.0f, -1.0f, -1.0f,  1.0f,  1.0f,
                     1.0f,  1.0f, 1.0f,  1.0f, -1.0f, -1.0f, -1.0f, -1.0f;
    cubeVertices *= (0.5f * hK_ * CUBE_SIZE);

    Eigen::MatrixXi cubeTriangles( 3, 12 );     // 12 triangles on the surface
    cubeTriangles << 0, 2, 1, 6, 7, 5, 4, 3, 4, 1, 3, 6,
                     1, 3, 5, 2, 6, 4, 0, 7, 5, 0, 2, 7,
                     2, 0, 6, 1, 5, 7, 3, 4, 1, 4, 6, 3;

    p.reserve( knots_p_.size() * 8 );
    t.reserve( knots_p_.size() * 12 );
    for (size_t i=0; i<signalingCenters_.size(); ++i)
    {
        Vec3f knot_p = signalingCenters_[i];
        // Add cube vertices, triangles.
        for (int j=0; j<8; ++j)
            p.push_back( cubeVertices.col(j) + knot_p );
        for (int j=0; j<12; ++j)
            t.push_back( cubeTriangles.col(j) + 8*i*Vec3i::Ones(3,1) );
    }
}



int Knots::getSourceNodes( const LevelSet& lset, 
                           std::vector<int>& sourceNodes ) const
{
    sourceNodes.clear();
    sourceNodes.reserve( signalingCenters_.size() );

    for (size_t i=0; i<signalingCenters_.size(); ++i) {
        Vec3f p = signalingCenters_[i];

        Vec3f normal_p;         // level set normal at p
        float lset_p;           // level set at p
        if (lset.interpolatePhi(p, lset_p, normal_p))
            return -1;

        // Interior source position.
        Vec3f src_p = p - normal_p * (lset_p + SOURCE_DISTANCE*hK_);

        // Search for nodes with hK_ radius of src_p.
        auto nodes = lset.getMesh().getNodesWithinRadius(src_p, hK_);
        if (nodes.size() == 0 || nodes[0] == -1) {
            std::cerr << __FUNCTION__ << ": Cannot determine source position "
                      << "for knot " << i << std::endl;
            return -1;
        }
        int node = nodes[0];    // nodes[0] is guaranteed to be closest to src_p.
        auto& phi = lset.getPhi();
        if (phi(node) > 0.0f) {
            std::cerr << "Error: Source node " << node << " outside the interface"
                      << " (phi = " << phi(node) << ")" << std::endl;
        }

        sourceNodes.push_back(node);
    }

    return 0;
}


//
// PRIVATE METHODS
//


int Knots::createSignalingCenters_( const LevelSet& lset,
                                    const std::vector<std::vector<int>>& clusters,
                                    const std::vector<Vec3f>& knots_p )
{
    for (auto& cluster : clusters) {
        Vec3f mean_p( 0.0f, 0.0f, 0.0f );
        for (const int& i : cluster)
            mean_p += knots_p[i];
        mean_p /= cluster.size();
       
        // Project the mean cluster position to the interface.
        Vec3f normal_p;     // level set normal at mean_p
        float lset_p;       // level set at mean_p
        if (lset.interpolatePhi(mean_p, lset_p, normal_p))
            return  -1;
        mean_p -= lset_p * normal_p;

        signalingCenters_.push_back( mean_p );
    }

    return 0;
}



void Knots::fillCluster_( const int cnode, const std::vector<Vec3f>& knots_p,
                          std::vector<int>& cluster, std::vector<int>& nodes, 
                          int depth ) const
{
    auto nodes_local = nodes;

    for (int node : nodes_local) {
        // Nodes that are within CLUSTER_DIAMETER * hK_ from each other will be 
        // taken as forming a single cluster.
        Vec3f d = knots_p[node] - knots_p[cnode];
        if (d.norm() < CLUSTER_DIAMETER * hK_) {
            cluster.push_back( node );
            nodes.erase( std::remove(nodes.begin(), nodes.end(), node),
                         nodes.end() );
        }
    }

    depth++;    // recursion depth
    if ((int)cluster.size() > depth)
        fillCluster_( cluster[depth], knots_p, cluster, nodes, depth );
}



void Knots::createKnotClusters_( const std::vector<Vec3f>& knots_p,
                                 std::vector<std::vector<int>>& clusters ) const
{
    // Fill node IDs
    std::vector<int> nodes( knots_p.size() );
    std::iota( nodes.begin(), nodes.end(), 0 );

    while (nodes.size() > 0) {
        int node = nodes[0];
        std::vector<int> cluster = {node};
        nodes.erase( std::remove(nodes.begin(), nodes.end(), node),
                     nodes.end() );

        fillCluster_( node, knots_p, cluster, nodes, 0 );
        clusters.push_back( cluster );
    }
}
