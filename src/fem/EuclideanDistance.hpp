//
// Resets the level function to a distance function by computing the minimum
// Euclidean distances from all mesh nodes to the interface.
//
// Uses CUDA device if available, otherwise runs on CPU.
//

#pragma once

#include <cuda_runtime.h>


class EuclideanDistance
{

public:
    EuclideanDistance() : 
        h_p(nullptr), 
        d_p(nullptr), 
        d_phi(nullptr), 
        m_meshSize(0),
        m_devCount(0),
        m_threads(1)
    {}

    //
    // Upload the fixed mesh node position (meshp, 3xn) and initial level 
    // set (phi0, n).
    //
    void        init( const float* meshp, const double* phi0, const int n,
                      const int nThreads );

    //
    // Given interface-mesh intersection points (ifp, size n_ifp), computes the 
    // signed distance function (phi_n, size m_meshSize).
    // Note that phi_n must contain the old (pre-reset) distances.
    //
    void        resetPhi( const float* ifp, const int n_ifp, double* phi_n );


private:

    const float*    h_p;        // host mesh node positions
    float*          d_p;        // GPU mesh node positions (x1, y1, z1, x2, ...)
    double*         d_phi;      // GPU level set function
    
    int             m_meshSize; // number of mesh nodes.

    int             m_devCount; // number of CUDA capable devices
    int             m_threads;  // max. number of threads for host-side operations

};
