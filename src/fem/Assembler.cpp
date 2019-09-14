#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Geometry>

#include "Assembler.hpp"
#include "MathTypes.hpp"
#include "common/Utils.hpp"


namespace FEM {

void Assembler::grad( const Mesh3D& mesh, const VecXd u, 
                      std::vector<Vec3f>& grad_t ) const
{
    grad_t.resize( mesh.t().size() );

    // 2nd order quadrature rule for P1 element:
    Quadrature qr(2, mesh.getElement(1));

    #pragma omp parallel num_threads( nThreads_ )
    {
    std::vector<Mat34d, Eigen::aligned_allocator<Mat34d>> gphi_global(4);

    #pragma omp for 
    for (size_t k=0; k<mesh.t().size(); ++k) {
        auto& t = mesh.t(k);

        double detB;
        Mat3d invBt;
        getElementInfo_ij_( mesh, k, detB, invBt );
        for (int j=0; j<4; j++)
            gphi_global[j] = invBt * qr.gphi()[j];
        // NOTE: Assuming here that gradients are constants over elements, i.e., 
        // col(i) = col(j), i,j € [0,3]. So we only compute the gradient at the 
        // first quadrature point (col(0)).
        Vec3d g = u(t(0)) * gphi_global[0].col(0) + 
                  u(t(1)) * gphi_global[1].col(0) +
                  u(t(2)) * gphi_global[2].col(0) + 
                  u(t(3)) * gphi_global[3].col(0);
        grad_t[k] = g.cast<float>();
    }
    }   // END pragma omp
}



void Assembler::grad_nodes( const Mesh3D& mesh, const VecXd u, 
                            std::vector<Vec3f>& grad_n ) const
{
    // Gradient interpolation is done by first finding all elements
    // associated with a node, then setting the node gradient as the average
    // of the element gradients weighted by element volumes.

    std::vector<Vec3f> grad_t;
    grad( mesh, u, grad_t );

    grad_n.resize( mesh.p().size() );

    // For each node compute the mean normal vector from element normals.
    #pragma omp parallel for num_threads( nThreads_ )
    for (size_t i=0; i<mesh.p().size(); ++i) {
        VecXi t_idx = mesh.p2t(i);      // indices of t's associated with node
        Vec3f n( 0.0f, 0.0f, 0.0f );    // mean normal

        // TODO: Currently assuming tetras are all equal size!
        for (int j=0; j<t_idx.size(); ++j)
            n += grad_t[ t_idx(j) ];
        n /= (float)t_idx.size();

        grad_n[i] = n;
    }
}



void Assembler::lump_mass( const sparse_matrix& M, sparse_matrix& invM ) const
{
    uint32_t n = M.rows();
    VecXd ones = VecXd::Ones(n,1);
    auto rsum = M * ones;
    VecXd d = ones.cwiseQuotient( rsum );
    sparse_matrix I(n,n);
    I.setIdentity();
    invM = d.asDiagonal() * I;      // much faster than calling sparseView()
}

}   // END namespace
