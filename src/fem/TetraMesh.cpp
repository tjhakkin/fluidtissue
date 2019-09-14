#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <cmath>
#include <chrono>

#include <eigen3/Eigen/Geometry>

#include <CGAL/Triangulation_3.h>
#include <tbb/task_scheduler_init.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include "MathTypes.hpp"
#include "Assembler.hpp"
#include "TetraMesh.hpp"


namespace FEM {

int Mesh3D::Delaunay_3D( const Parameters& parameters, const int nThreads )
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
    typedef CGAL::Triangulation_vertex_base_with_info_3< unsigned, K >  Vb;
    typedef CGAL::Triangulation_cell_base_3< K >                        Cb;
    typedef CGAL::Triangulation_data_structure_3< Vb, Cb, CGAL::Parallel_tag >  Tds;
    typedef CGAL::Spatial_lock_grid_3< CGAL::Tag_priority_blocking >            Lock_ds;
    typedef CGAL::Delaunay_triangulation_3< K, Tds, CGAL::Default, Lock_ds >    Delaunay;

    // Set up vertices for triangulation.
    std::vector< std::pair<Delaunay::Point, unsigned> > points( p_.size() );
    for (uint32_t i=0; i<p_.size(); i++) {
        Vec3d v = p_.at(i).cast<double>();
        points[i] = { Delaunay::Point(v(0), v(1), v(2)), i };
    }

    std::cout << "  - Triangulating... " << std::flush;
    tbb::task_scheduler_init init( nThreads );
    // Unit bounding box for lock grid:
    // CGAL::Bbox_3 bbox(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
    // TODO: N_max should be "cells per axis" per the CGAL documentation;
    // presumably doesn't need to be the exact value, but what are the limits?
    int density = (int)parameters["Node density"];
    int N_max = 2 * std::floor( std::cbrt(density + 1) );
    Lock_ds locking_ds( CGAL::Bbox_3(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0), N_max );

    Delaunay dt;
    dt.set_lock_data_structure( &locking_ds );
    dt.insert( points.begin(), points.end() );
    // Delaunay dt( points.begin(), points.end(), &locking_ds );
    // Delaunay dt( points.begin(), points.end() );         // single-threaded
    std::cout << "done." << std::endl;

    //
    // Fill mesh.edges
    //
    std::map< std::pair<int,int>, int > edgeIndices;
    int n_edges = std::distance( dt.finite_edges_begin(), dt.finite_edges_end() );
    edges_.reserve( n_edges );

    std::cout << "  - Filling mesh.edges... " << std::flush;
    for (Delaunay::Finite_edges_iterator it = dt.finite_edges_begin();
         it != dt.finite_edges_end(); ++it)
    {
        Vec2i e = { (int)it->first->vertex( (it->second) )->info(),
                    (int)it->first->vertex( (it->third) )->info() };
        edges_.push_back( e );          // TODO: Add sorted edges?

        if (e(0) < e(1))
            edgeIndices[ {e(0),e(1)} ] = edgeIndices.size();
        else
            edgeIndices[ {e(1),e(0)} ] = edgeIndices.size();
    }
    std::cout << "done." << std::endl;

    //
    // Fill mesh.t (tetrahedra nodes), mesh.t2e (tetra-edges maps)
    //
    int n_cells = std::distance( dt.finite_cells_begin(), dt.finite_cells_end() );
    t_.reserve( n_cells );
    t2e_.reserve( n_cells );
    // For storing all node-tetrahedron pairs; needed later for filling mesh.p2t:
    std::vector< std::pair<int,int> > pt;
    pt.reserve(4 * n_cells);

    std::cout << "  - Filling mesh.t, mesh.t2e... " << std::flush;
    for (Delaunay::Finite_cells_iterator it = dt.finite_cells_begin();
         it != dt.finite_cells_end(); ++it)
    {
        Vec4i t = { (int)it->vertex(0)->info(), (int)it->vertex(1)->info(),
                    (int)it->vertex(2)->info(), (int)it->vertex(3)->info() };
        t_.push_back(t);

        // node-tetrahendon pairs (needed later)
        int ti = t_.size() - 1;
        for (int j=0; j<t.size(); j++)
            pt.push_back( {t(j), ti} );

        Vec6i t2e;
        // 6 edges per tetra, node index permutations:
        std::vector<Vec2i> P( {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}} );
        for (int j=0; j<P.size(); j++) {
            std::pair<int,int> edge = { t(P[j](0)), t(P[j](1)) };
            if (edge.second < edge.first)
                edge = { edge.second, edge.first };
            t2e(j) = edgeIndices[edge];
        }
        t2e_.push_back( t2e );
    }
    std::cout << "done." << std::endl;

    //
    // Fill mesh.p2t (node-to-elements maps)
    //
    std::sort(pt.begin(), pt.end());        // sorts according to pt.first
    p2t_.resize( p_.size() );

    std::cout << "  - Filling mesh.p2t... " << std::flush;
    int start = 0;
    for (size_t i=1; i<pt.size()+1; i++) {
        int node = pt[start].first;
        if (i == pt.size() || pt[i].first != node) {
            VecXi p2t(i - start);
            for (int j=start; j<i; j++)
                p2t(j - start) = pt[j].second;
            p2t_[node] = p2t;
            start = i;
        }
    }
    std::cout << "done." << std::endl;

    //
    // Fill mesh.b_tris (outer/exterior boundary facets, oriented), 
    // mesh.b_edgeIDs (boundary facet edge indices)
    //
    std::vector<Delaunay::Cell_handle> facets;
    dt.incident_cells( dt.infinite_vertex(), std::back_inserter(facets) );
    b_tris_.reserve( facets.size() ); 
    b_edgeIDs_.reserve( facets.size() );
    
    std::cout << "  - Filling mesh.b_tris, mesh.b_edgeIDs... " << std::flush;
    for (auto f : facets) {
        int inf_vert = f->index( dt.infinite_vertex() );
        bool even = inf_vert%2==0;
        Vec3i t = { (int)f->vertex( (inf_vert + (even ? 3:1))%4 )->info(),
                    (int)f->vertex( (inf_vert + 2)%4 )->info(),
                    (int)f->vertex( (inf_vert + (even ? 1:3))%4 )->info() };
        b_tris_.push_back(t);

        Vec3i edgeIDs;
        for (int j=0; j<edgeIDs.size(); j++) {
            std::pair<int,int> edge = { t(j%3), t((j+1)%3) };   // indices 01, 12, 20
            if (edge.second < edge.first)
                edge = { edge.second, edge.first };
            edgeIDs(j) = edgeIndices[edge];
        }
        b_edgeIDs_.push_back( edgeIDs );
    }

    std::cout << "done." << std::endl;

    return 0;
}



void Mesh3D::initMesh( const Parameters& parameters, Vec3i& nodes, Vec3f& dimensions )
{
    // The final domain width, height, depth have the upper bounds (see below)
    // twice the given values, e.g., width = 1.0 indicates domain in range [-1,1] 
    // along the x axis.
    // Density is the total number of nodes within half-closed cube [0,1)^3.
    double width = parameters["Domain width"];
    double height = parameters["Domain height"];
    double depth = parameters["Domain depth"];
    int density = (int)parameters["Node density"];

    // Numbers of nodes along each axis (width, height, depth).
    int ns = std::floor( std::cbrt(density + 1) );  // +1 to avoid rounding under
    int nw = std::floor( 2.0*width*(ns+1) );
    int nh = std::floor( 2.0*height*(ns+1) );
    int nd = std::floor( 2.0*depth*(ns+1) );

    // Final, origin-centered domain dimensions.

    float h = 2.0 * std::max({width, height, depth}) / (std::max({nw, nh, nd}) - 1);
    width = 0.5f * (nw-1) * h;
    height = 0.5f * (nh-1) * h;
    depth = 0.5f * (nd-1) * h;

    // Output
    nodes = {nw, nh, nd};
    dimensions = {width, height, depth};

    // Generated mesh nodes.
    p_.reserve( nodes(0) * nodes(1) * nodes(2) );
    for (int i=0; i<nodes(0); i++) {
        float px = -dimensions(0) + i*h;
        // std::cerr << px << std::endl;

        for (int j=0; j<nodes(1); j++) {
            float py = -dimensions(1) + j*h;

            for (int k=0; k<nodes(2); k++) {
                float pz = -dimensions(2) + k*h;
                p_.push_back( Vec3f(px, py, pz) );
            }
        }
    }
}



void Mesh3D::initDofs()
{
    //
    // Fill tetrahedron -> dofs map.
    //
    VecXi dofI(element_[1].N().size());
    dofIndices_.resize(t_.size());

    for (size_t k=0; k<t_.size(); k++) {
        for (int i=0; i<t_[k].size(); i++)
            dofI(i) = t_[k](i);

        // Edge node-pair permutations. The order depends on the order in
        // which the basis functions are listed.
        // This could be simplified if e2t could be assumed to follow the same
        // order as the shape functions - assuming now that's not the case.

        // element edges as node pairs, as pulled from t2e
        std::vector<Vec2i> tedges(t2e_[k].size());
        for (size_t i=0; i<tedges.size(); i++) {
            int j = t2e_[k](i);
            tedges[i] = {edges_[j](0), edges_[j](1)};
            if (tedges[i](1) < tedges[i](0))
                tedges[i] = {edges_[j](1), edges_[j](0)};
        }

        // get the global element k dof edges using the reference element
        // dof indexing; find the corresponding global edge index for dofI.
        auto& en = element_[1].getEdgeDofs(); // reference element dof edges
        for (size_t i=0; i<en.size(); i++) {
            auto& e = en[i];
            Vec2i dofe = {t_[k](e(0)), t_[k](e(1))};
            if (dofe(1) < dofe(0))
                dofe = {t_[k](e(1)), t_[k](e(0))};

            // find the indices of the edges corresponding to the dof edges
            int j = std::find(tedges.begin(), tedges.end(), dofe) - tedges.begin();
            // if (j >= t2e[k].size())  // TODO: add error
            int n = t_[k].size();
            dofI(n+i) = t2e_[k](j) + p_.size();
        }

        dofIndices_[k] = dofI;
    }

    //
    // Fill list of interior dofs
    //

    // obtain unique node dofs
    std::vector<int> IB(p_.size());      // all node dofs (interior + boundary)
    std::iota(IB.begin(), IB.end(), 0);
    std::vector<int>nodeDofs = unique_nodes(b_tris_);
    interiorDofs_[0] = set_diff(IB, nodeDofs);  // interior node dofs

    // obtain unique edge dofs
    std::vector<int> edgeDofs;
    edgeDofs.reserve(3 * b_edgeIDs_.size());
    for (auto& te : b_edgeIDs_) {
        for (int i=0; i<3; i++) {
            // edge dof indexing begins at p.size():
            edgeDofs.push_back(te(i) + p_.size());
        }
    }
    auto last = std::unique(edgeDofs.begin(), edgeDofs.end());
    edgeDofs.erase(last, edgeDofs.end());

    size_t n = nodeDofs.size() + edgeDofs.size();
    std::vector<int>B(n);               // boundary dofs
    std::copy(nodeDofs.begin(), nodeDofs.end(), B.begin());
    std::copy(edgeDofs.begin(), edgeDofs.end(), B.begin()+nodeDofs.size());

    IB.resize(getNDofs(2));             // all dofs
    std::iota(IB.begin(), IB.end(), 0);

    interiorDofs_[1] = set_diff(IB, B);    // interior dofs
}



int Mesh3D::createRegularMesh( const Parameters& parameters, int nThreads )
{
    Vec3i nodes;        // number of nodes along each axis
    Vec3f dimensions;   // origin-centered domain dimension (width, height, depth)
    initMesh( parameters, nodes, dimensions );

    std::cout << "  - Mesh: " << nodes(0) << " x " << nodes(1) << " x " << nodes(2) << " within "
              << "[" << -dimensions(0) << "," << dimensions(0) << "] x "
              << "[" << -dimensions(1) << "," << dimensions(1) << "] x "
              << "[" << -dimensions(2) << "," << dimensions(2) << "]" << std::endl;

    Delaunay_3D( parameters, nThreads );

    initDofs();
    std::cout << "    Number of elements: " << t_.size() << std::endl
              << "    Number of nodes: " << p_.size() << std::endl;

    return 0;
}

}   // END namespace
