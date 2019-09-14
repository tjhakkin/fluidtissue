#include <cfloat>
#include <iostream>
#include "EuclideanDistance.hpp"


namespace {

static const int THREADS_PER_BLOCK = 512;


//
// Computes the shortest distance from mesh node (p; 3) to the interface (ifp; 3xn).
//
__host__ __device__
float distanceToInterface_( const float* p, const float* ifp, const int n_ifp )
{
    float dmin = FLT_MAX;
    
    for (size_t j=0; j<n_ifp; j++) {
        float dx = ifp[3*j+0] - p[0];
        float dy = ifp[3*j+1] - p[1];
        float dz = ifp[3*j+2] - p[2];

#ifdef __CUDA_ARCH__
        float d = norm3df( dx, dy, dz );
#else
        // Note: std::hypot() would be better, but requires C++17.
        float d = std::sqrt( dx*dx + dy*dy + dz*dz );
#endif
        if (d < dmin)
            dmin = d;
    }
    
    return dmin;
}

}   // END namespace


//
// Kernels begin
//


__global__ 
void distanceToInterface_k( const float* meshp, const int n_meshp, 
                            const float* ifp, int n_ifp, double* phi )
{
    int k = threadIdx.x + blockIdx.x*blockDim.x;
    if (k >= n_meshp)
        return;

    float p[3] = {meshp[3*k], meshp[3*k+1], meshp[3*k+2]};
    double dmin = distanceToInterface_( p, ifp, n_ifp );
    if (phi[k] < 0.0f)
        dmin = -dmin;
    phi[k] = dmin;
}


//
// Kernels end
//


void EuclideanDistance::init( const float* meshp, const double* phi0, const int n,
                              const int nThreads )
{
    m_meshSize = n;
    m_threads = nThreads;

    cudaGetDeviceCount( &m_devCount );

    if (m_devCount) {
        int size = 3 * n * sizeof(float);
        cudaMalloc( &d_p, size );
        cudaMemcpy( d_p, meshp, size, cudaMemcpyHostToDevice );

        size = n * sizeof(double);
        cudaMalloc( &d_phi, size );
        cudaMemcpy( d_phi, phi0, size, cudaMemcpyHostToDevice );
    }
    else {
        h_p = meshp;
    }
}



void EuclideanDistance::resetPhi( const float* ifp, const int n_ifp, double* phi_n )
{
    if (m_devCount) {
        int blockSize = THREADS_PER_BLOCK;
        int gridSize = (int)ceil( (float)m_meshSize / blockSize );

        float* d_ifp;
        int size = 3 * n_ifp * sizeof(float);
        cudaMalloc( &d_ifp, size );
        cudaMemcpy( d_ifp, ifp, size, cudaMemcpyHostToDevice );

        size = m_meshSize * sizeof(double);
        cudaMemcpy( d_phi, phi_n, size, cudaMemcpyHostToDevice );
        distanceToInterface_k<<< gridSize, blockSize >>>( d_p, m_meshSize, 
                                                          d_ifp, n_ifp, d_phi );   
        cudaFree( d_ifp );
        cudaMemcpy( phi_n, d_phi, size, cudaMemcpyDeviceToHost );
    }
    else {
        #pragma omp parallel for num_threads( m_threads )
        for (size_t i=0; i<m_meshSize; i++) {
            float p[3] = {h_p[3*i], h_p[3*i+1], h_p[3*i+2]};
            float dmin = distanceToInterface_( p, ifp, n_ifp );
            if (phi_n[i] < 0.0f)
                dmin = -dmin;
            phi_n[i] = dmin;
        }
    }
}
