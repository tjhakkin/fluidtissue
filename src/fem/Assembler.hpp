//
// FEM matrix assembler
//
// ** Basic usage **
// 
// Initialize the assembler and assemble the standard mass matrix (u.*v) with
//
// FEM::Assembler fem;
// auto wM = WEAK_FORM( u.cwiseProduct(v) )
// FEM::sparse_matrix M; 
// fem.assemble_bilin( mesh, wM, {}, M );
//
// Assemble the quadratic (ux.*vx) with per-element data 
// 
// VecXd data = ...;
// auto wA11 = WEAK_FORM( ux.cwiseProduct(vx) * Ve )
// fem.assemble_bilin( mesh, wA11, data, A11, 2, 2 );
//
// where 'Ve' is the keyword for the assembler to distribute the values of data
// per element.
//
// TODO:
// assemble_bilin() currently calls either bilin_without_p2t(), bilin_with_p2t() 
// or bilin_without_p2t_ij(), depending on the situation. The _ij variant is more 
// general and should be the only one used, but need to implement support for 
// looping over nodes instead of elements and support for 3D data first.
//

#pragma once

#include <thread>

#include "TetraMesh.hpp"
#include "Tetrahedron.hpp"
#include "Quadrature.hpp"
#include "MathTypes.hpp"


// Constructs weak formulation lambda for assemble_bilin().
#define WEAK_FORM( handle )                                                                                                         \
[&](const Eigen::Ref<const VecXd>& u, const Eigen::Ref<const VecXd>& v,                                                             \
    const Eigen::Ref<const VecXd>& vx, const Eigen::Ref<const VecXd>& vy, const Eigen::Ref<const VecXd>& vz,                        \
    const Eigen::Ref<const VecXd>& ux, const Eigen::Ref<const VecXd>& uy, const Eigen::Ref<const VecXd>& uz,                        \
    float h, double Ve, const Eigen::Ref<const VecXd>& V1, const Eigen::Ref<const VecXd>& V2, const Eigen::Ref<const VecXd>& V3)    \
{                                                                                                                                   \
    return handle;                                                                                                                  \
};


namespace FEM {

class Assembler {

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Assembler() : 
        nThreads_(std::thread::hardware_concurrency())
    {}
    
    //
    // Assembles a general bilinear form A for given Lambda f for polynomial 
    // degrees i, j, with optional element data (e.g., viscosities).
    //
    template<typename Lambda>
    int assemble_bilin( const FEM::Mesh3D& mesh, Lambda& f, const VecXd& data, 
                        sparse_matrix& A, int degree_i = 1, int degree_j = 1 ) const
    {
        size_t nt = mesh.t().size();
        size_t np = mesh.p().size();
        if (data.size() != 0 && data.size() != nt && data.size() != 3*np) {
            std::cerr << __FUNCTION__ << "(): Incorrect data vector size " << data.size()
                      << ". Must be either 0, number of elements (" << nt 
                      << ") or 3 x number of nodes (" << 3*np << ")." << std::endl;
            return -1;
        }

        // If mesh contains p2t (nodes-to-tetrahedra), loops over nodes and 
        // fill the sparse matrix row at a time. Doesn't require much memory.
        // Support only P1-P1 element.
        if (mesh.p2t().size() == np && degree_i == 1 && degree_j == 1) {
            // Average number of non-zeros per row of A = edges / nodes * 2 + 1,
            // where x2 since two nodes per edge and +1 to add the row node (dianonal).
            int nnz = std::ceil( 2.0 * mesh.edges().size() / mesh.p().size() + 1 );
            A.reserve( Eigen::VectorXi::Constant(mesh.p().size(), nnz) );
            bilin_with_p2t( mesh, f, data, A );
        }
        else if (degree_i == 1 && degree_j == 1) {
            bilin_without_p2t( mesh, f, data, A );
        }
        else if (degree_i >= 1 && degree_j >= 1) {
            // Loops over elements; uses more memory and doesn't support 3D data.
            // Supports both P1 and P2 elements.
            // Intended to replace the other two assemblers.
            bilin_without_p2t_ij( mesh, degree_i, degree_j, f, data, A );
        }
        else {
            std::cerr << __FUNCTION__ << "(): Invalid polynomial degrees given, "
                      << "both must be >= 1." << std::endl;
            return -1;
        }

        return 0;
    };

    //
    // Compute element gradients of given field u.
    //
    void grad( const FEM::Mesh3D& mesh, const VecXd u, 
               std::vector<Vec3f>& grad_t ) const;

    //
    // Compute gradients of given field u at mesh nodes.
    //
    // Limitations:
    // - Assumes all tetrahedra are equal size.
    //
    void grad_nodes( const Mesh3D& mesh, const VecXd u, 
                     std::vector<Vec3f>& grad_n ) const;

    //
    // Approximates the inverse of mass by performing mass lumping.
    //
    void lump_mass( const sparse_matrix& M, sparse_matrix& invM ) const;


private:
    
    //
    // Called from assemble_bilin() to assemble the bilinear form P1-P1 for
    // function handle f with p2t information in mesh.
    //
    template<typename Lambda>
    void bilin_with_p2t( const FEM::Mesh3D& mesh, Lambda& f, const VecXd& data, 
                         sparse_matrix& A ) const
    {
        size_t nt = mesh.t().size();
        size_t np = mesh.p().size();

        // 2nd order quadrature rule for P1 element:
        Quadrature qr(2, mesh.getElement(1));
        // quadrature weights, scale by tetrahedra volume 1/6:
        VecXd qw = qr.qw() / 6.0;

        #pragma omp parallel num_threads( nThreads_ )
        {
        std::vector<Vec4d, Eigen::aligned_allocator<Vec4d>> uv_x(4);
        std::vector<Vec4d, Eigen::aligned_allocator<Vec4d>> uv_y(4);
        std::vector<Vec4d, Eigen::aligned_allocator<Vec4d>> uv_z(4);

        #pragma omp for
        for (size_t i=0; i<np; i++) {
            std::map<int, double> map;      // column, integral value

            for (size_t ti=0; ti<mesh.p2t(i).size(); ti++) {
                int k = mesh.p2t(i)(ti);    // element containing node i
                auto& t = mesh.t(k);        // element k nodes

                // Compute affine matrix determinant, basis gradients.
                double detB;
                Mat3d invBt;
                getElementInfo_ij_( mesh, k, detB, invBt );
                for (int j=0; j<4; j++) {
                    Mat34d gphi_global = invBt * qr.gphi()[j];
                    uv_x[j] = gphi_global.row(0);
                    uv_y[j] = gphi_global.row(1);
                    uv_z[j] = gphi_global.row(2);
                }

                // Per-element data (1D, e.g., viscosity)
                double Ve = 1.0;
                if (data.size() == nt)
                    Ve = data[k];

                // Per-node data (3D, e.g., velocity field)
                Vec4d V1(1.0, 1.0, 1.0, 1.0);
                Vec4d V2(1.0, 1.0, 1.0, 1.0);
                Vec4d V3(1.0, 1.0, 1.0, 1.0);
                if (data.size() == 3*np) {
                    // Computing the contributions of the node values at
                    // quadrature points.
                    V1.setZero();
                    V2.setZero();
                    V3.setZero();
                    for (int j=0; j<4; j++) {
                        V1 += data[t(j)+0*np] * qr.phi()[j];
                        V2 += data[t(j)+1*np] * qr.phi()[j];
                        V3 += data[t(j)+2*np] * qr.phi()[j];
                    }
                }

                // Mesh parameter, proportional to the element diameter:
                float h = std::cbrt( std::fabs(detB) );

                // Evaluate v, v' at the quadrature points; first need to
                // find node i in the element k:
                int j = 0;
                if (t(1) == i)
                    j = 1;
                else if (t(2) == i)
                    j = 2;
                else if (t(3) == i)
                    j = 3;
                else {}
                Vec4d v = qr.phi()[j];
                Vec4d vx = uv_x[j];
                Vec4d vy = uv_y[j];
                Vec4d vz = uv_z[j];

                // Integrate over element k.
                for (int j=0; j<4; j++) {
                    Vec4d u = qr.phi()[j];
                    Vec4d ux = uv_x[j];
                    Vec4d uy = uv_y[j];
                    Vec4d uz = uv_z[j];

                    Vec4d feval = f(u,v,vx,vy,vz,ux,uy,uz,h,Ve,V1,V2,V3);
                    double qint = feval.dot(qw) * detB;
                    map[t(j)] += qint;
                }
            }

            for (auto& j : map)
                A.insert(i, j.first) = j.second;
        }
        }   // END pragma omp

        A.makeCompressed();
    }

    //
    // Called from assemble_bilin() to assemble the bilinear form P1-P1 for
    // function handle f without p2t information in mesh.
    //
    // TODO: Implement support for 3D data (e.g., for advection velocities).
    //
    template<typename Lambda>
    void bilin_without_p2t( const FEM::Mesh3D& mesh, Lambda& f, const VecXd& data, 
                            sparse_matrix& A ) const
    {
        size_t nt = mesh.t().size();

        typedef Eigen::Triplet<double> T;
        std::vector<T> triplets_A;
        // 4 basis functions => 4^2 = 16 pairs per element in the integration:
        triplets_A.resize(16*nt);

        // 2nd order quadrature rule for P1 element:
        Quadrature qr(2, mesh.getElement(1));
        // quadrature weights, scale by tetrahedra volume 1/6:
        VecXd qw = qr.qw() / 6.0;

        #pragma omp parallel
        {
        std::vector<Vec4d, Eigen::aligned_allocator<Vec4d>> uv_x(4);
        std::vector<Vec4d, Eigen::aligned_allocator<Vec4d>> uv_y(4);
        std::vector<Vec4d, Eigen::aligned_allocator<Vec4d>> uv_z(4);

        #pragma omp for 
        for (uint32_t k=0; k<nt; k++) {
            auto& t = mesh.t(k);

            // Compute affine matrix determinant, basis gradients.
            double detB;
            Mat3d invBt;
            getElementInfo_ij_( mesh, k, detB, invBt );
            for (int j=0; j<4; j++) {
                Mat34d gphi_global = invBt * qr.gphi()[j];
                uv_x[j] = gphi_global.row(0);
                uv_y[j] = gphi_global.row(1);
                uv_z[j] = gphi_global.row(2);
            }
            
            // Per-element data (1D, e.g., viscosity)
            double Ve = 1.0;
            if (data.size() == nt)
                Ve = data[k];

            // Per-node data (3D, e.g., velocity field)
            Vec4d V1(1.0, 1.0, 1.0, 1.0);
            Vec4d V2(1.0, 1.0, 1.0, 1.0);
            Vec4d V3(1.0, 1.0, 1.0, 1.0);
            // TODO: implement

            // Mesh parameter, proportional to the element diameter:
            float h = std::cbrt( std::fabs(detB) );

            // Integrate over element k.
            for (int i=0; i<4; i++) {
                Vec4d u = qr.phi()[i];
                Vec4d ux = uv_x[i];
                Vec4d uy = uv_y[i];
                Vec4d uz = uv_z[i];

                for (int j=0; j<4; j++) {
                    Vec4d v = qr.phi()[j];
                    Vec4d vx = uv_x[j];
                    Vec4d vy = uv_y[j];
                    Vec4d vz = uv_z[j];

                    Vec4d feval = f(u,v,vx,vy,vz,ux,uy,uz,h,Ve,V1,V2,V3);
                    double qint = feval.dot(qw) * detB;
                    uint32_t n = 16*k + i*4 + j;
                    triplets_A[n] = T(t(j), t(i), qint);
                }
            }
        }
        }   // END pragma omp

        A.setFromTriplets( triplets_A.begin(), triplets_A.end() );
    }

    //
    // Called from assemble_bilin() to assemble the bilinear form Pi-Pj for
    // given polynomial degrees i, j, with function handle f and without p2t
    // information in mesh.
    //
    // TODO: Implement support for 3D data (e.g., for advection velocities).
    //    
    template<typename Lambda>
    void bilin_without_p2t_ij( const FEM::Mesh3D& mesh, int degree_i, int degree_j, 
                               Lambda& f, const VecXd& data, sparse_matrix& A ) const
    {
        // Required quadrature degree = degree_i + degree_j. For computational,
        // efficiency, always use the lowest order method possible.

        // obtain quadrature rules 
        int d = degree_i + degree_j;
        Quadrature qr_i(d, mesh.getElement(degree_i));
        Quadrature qr_j(d, mesh.getElement(degree_j));

        // quadrature weights, scale by tetrahedra volume 1/6.
        // NOTE: qr_i, qr_j were initialized with the same degree.
        VecXd qw = qr_i.qw() / 6.0;

        // reference elements gradients:
        auto& gphi_i = qr_i.gphi();
        auto& gphi_j = qr_j.gphi();
        int nDofI = gphi_i.size();
        int nDofJ = gphi_j.size();
        size_t nt = mesh.t().size();

        int ni = mesh.getNDofs(degree_i);
        int nj = mesh.getNDofs(degree_j);
        A = sparse_matrix(ni, nj);

        // If memory consumption is an issue, use a value of jobSize < nt;
        // splits the assembly into according sized pieces.
        int jobSize = nt;

        typedef Eigen::Triplet<double> T;
        std::vector<T> triplets_A;
        triplets_A.resize(nDofI * nDofJ * jobSize); // total number of integrations

        #pragma omp parallel num_threads( nThreads_ )
        {
        std::vector<VecXd> uv_x(nDofJ + nDofI);
        std::vector<VecXd> uv_y(nDofJ + nDofI);
        std::vector<VecXd> uv_z(nDofJ + nDofI);

        #pragma omp for
        for (uint32_t k=0; k<nt; k++) {
            // first 4 dofs are degree 1, all 10 degree 2:
            VecXi dofs = mesh.getDofIndices(k); 

            // Compute affine matrix determinant, basis gradients.
            double detB;
            Mat3d invBt;
            getElementInfo_ij_( mesh, k, detB, invBt );
            
            // compute gradient components at quadrature points
            for (int j=0; j<nDofJ; j++) {
                MatX gphi = invBt * gphi_j[j];
                uv_x[j] = gphi.row(0);
                uv_y[j] = gphi.row(1);
                uv_z[j] = gphi.row(2);
            }
            for (int i=0; i<nDofI; i++) {
                MatX gphi = invBt * gphi_i[i];
                uv_x[nDofJ + i] = gphi.row(0);
                uv_y[nDofJ + i] = gphi.row(1);
                uv_z[nDofJ + i] = gphi.row(2);
            }

            // Per-element data (1D, e.g., viscosity)
            double Ve = 1.0;
            if (data.size() == nt)
                Ve = data[k];

            // Per-node data (3D, e.g., velocity field)
            VecXd V1 = VecXd::Ones(qw.size());
            VecXd V2 = VecXd::Ones(qw.size());
            VecXd V3 = VecXd::Ones(qw.size());
            // TODO: implement

            // Mesh parameter, proportional to the element diameter:
            float h = std::cbrt( std::fabs(detB) );

            // Integrate over element k.
            for (int j=0; j<nDofJ; j++) {
                const VecXd& u = qr_j.phi()[j];
                VecXd& ux = uv_x[j];
                VecXd& uy = uv_y[j];
                VecXd& uz = uv_z[j];

                for (int i=0; i<nDofI; i++) {
                    const VecXd& v = qr_i.phi()[i];
                    VecXd& vx = uv_x[nDofJ + i];
                    VecXd& vy = uv_y[nDofJ + i];
                    VecXd& vz = uv_z[nDofJ + i];

                    VecXd feval = f(u,v,vx,vy,vz,ux,uy,uz,h,Ve,V1,V2,V3);
                    double qint = feval.dot(qw) * detB;

                    uint32_t n = nDofI*nDofJ*(k) + j*nDofI + i;
                    triplets_A[n] = T(dofs(i), dofs(j), qint);
                }
            }
        }
        }   // END pragma omp
        
        A.setFromTriplets( triplets_A.begin(), triplets_A.end() );

    }

    //
    // Given mesh and element index k, computes determinant of affine matrix A
    // and basis function gradients for k.
    //
    inline int getElementInfo_ij_( const FEM::Mesh3D& mesh, const int k, 
                                   double& detB, Mat3d& invBt ) const
    {
        auto& t = mesh.t(k);
        Vec3d p1 = mesh.p(t(0)).cast<double>();
        Vec3d p2 = mesh.p(t(1)).cast<double>();
        Vec3d p3 = mesh.p(t(2)).cast<double>();
        Vec3d p4 = mesh.p(t(3)).cast<double>();

        // Affine map Bx + b from local to global element; just need A:
        Mat3d B(3,3);
        B.col(0) = p2 - p1;
        B.col(1) = p3 - p1;
        B.col(2) = p4 - p1;

        // Determinant and inverse transpose.
        detB = std::fabs( B.determinant() );
        if (detB == 0) {
            std::cerr << "Singular matrix in " << __FUNCTION__ << ":"
                      << std::endl;
            std::cerr << std::endl << "B: " << std::endl << B << std::endl;
            std::cerr << std::endl << "Element indices: " << t(0) << ", "
                      << t(1) << ", " << t(2) << ", " << t(3) << std::endl;
            std::cerr << std::endl << "p1: " << p1.transpose() << std::endl
                                   << "p2: " << p2.transpose() << std::endl
                                   << "p3: " << p3.transpose() << std::endl;
            return -1;
        }

        invBt = B.inverse().transpose();

        return 0;
    }


    int     nThreads_;  // max. number of CPU threads
};

}   // END namespace
