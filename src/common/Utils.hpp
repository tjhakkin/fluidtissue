//
// Common helper functions for data sorting, timing etc.
//

#pragma once

#include <iostream>
#include <algorithm>
#include <chrono>


namespace FEM {

typedef std::chrono::steady_clock::time_point       chrono_time;

//
// Sorts given data vector v, returns sort indices.
//
template <typename T>
std::vector<int> sort_indices( T& v )
{
    // Sort indices
    std::vector<int> idx( v.size() );
    std::iota( idx.begin(), idx.end(), 0 );
    std::sort( idx.begin(), idx.end(),
               [&](int i, int j) { return v[i] < v[j]; } );

    // Sort data
    std::sort( v.begin(), v.end() );

    return idx;
}

//
// Returns set difference of two vectors.
//
template <typename T>
std::vector<T> set_diff( std::vector<T> set1, std::vector<T> set2 )
{
    std::sort( set1.begin(), set1.end() );
    std::sort( set2.begin(), set2.end() );

    std::vector<T> diff;
    std::set_difference( set1.begin(), set1.end(), set2.begin(), set2.end(),
                         std::inserter(diff, diff.begin()) );

    return diff;
}

//
// Returns the list of unique nodes found in elements (vector of Eigen vectors). 
//
template <typename T>
std::vector<int> unique_nodes( const std::vector<T>& elements )
{
    std::vector<int> nodes(0);
    if (elements.size() == 0)
        return nodes;

    int len = elements[0].size();
    nodes.reserve( elements.size() * len );

    for (auto& e : elements)
        nodes.insert( nodes.end(), e.data(), e.data()+len );
    std::sort( nodes.begin(), nodes.end() );
    auto last = std::unique( nodes.begin(), nodes.end() );
    nodes.erase( last, nodes.end() );

    return nodes;
}

//
// Record current time.
//
inline void tic( chrono_time& now )
{
    now = std::chrono::steady_clock::now();
}

//
// Returns time elapsed since 'start'.
//
inline double toc( chrono_time& start )
{
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    auto val = std::chrono::duration<double,std::milli>(diff).count();
    return val;
}

}   // END namespace
