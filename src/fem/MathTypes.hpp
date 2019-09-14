//
// Defines shorthard for various vector and matrix types provided by the Eigen
// library. Simple sparse submatrix implemention included.
//

#pragma once

#include <iostream>
#include <numeric>
#include <vector>

#include <eigen3/Eigen/Sparse>


namespace FEM {

typedef Eigen::Matrix<double, 3, Eigen::Dynamic>         MatX;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>   matrix;
// RowMajor sparse matrix seems to be a bit faster for solving RD at least:
typedef Eigen::SparseMatrix<double, Eigen::RowMajor>            sparse_matrix;
typedef std::vector<matrix>                                     matrix_array;

// Care needed with types that are multiples of 16 bytes to avoid alignment issues.
// Problematic types marked below. For more info, see
// https://eigen.tuxfamily.org/dox/group__TopicFixedSizeVectorizable.html
typedef Eigen::VectorXd                     VecXd;
typedef Eigen::VectorXf                     VecXf;
typedef Eigen::VectorXi                     VecXi;

typedef Eigen::Vector3d                     Vec3d;
typedef Eigen::Vector4d                     Vec4d;  // problem 
typedef Eigen::Vector2f                     Vec2f;
typedef Eigen::Vector3f                     Vec3f;
typedef Eigen::Vector4f                     Vec4f;  // problem
typedef Eigen::Vector2i                     Vec2i;
typedef Eigen::Vector3i                     Vec3i;
typedef Eigen::Vector4i                     Vec4i;  // problem
typedef Eigen::Matrix<int, 6, 1>            Vec6i; 
typedef Eigen::Matrix<double, 32, 1>        Vec32d; // problem
typedef Eigen::Matrix<double, 16, 1>        Vec16d; // problem

typedef Eigen::Matrix<float, 3,3>           Mat3f;
typedef Eigen::Matrix<double, 3,3>          Mat3d;
typedef Eigen::Matrix<float, 3,4>           Mat34f; // problem
typedef Eigen::Matrix<double, 3,4>          Mat34d; // problem


//
// Returns submatrix A(rows,cols) for row major matrix A.
//
inline void subSparse( sparse_matrix& A, const std::vector<int>& rows, 
                       const std::vector<int>& cols )
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletsB;
    tripletsB.reserve(A.nonZeros());

    // Submatrix column indices.
    std::vector<int> subCol(A.cols(), -1);
    for (size_t i=0; i<cols.size(); i++)
        subCol[cols[i]] = i;

    int row = 0;
    for (int i : rows) {
        // Iterate over columns of row i (rowMajor).
        // Note the order in which the columns are iterated is not known. 
        for (sparse_matrix::InnerIterator it(A,i); it; ++it) {
            int col = subCol[it.col()];
            if (col < 0)
                continue;
            tripletsB.push_back( T(row, col, it.value()) );
        }
        row++;
    }

    A = sparse_matrix(rows.size(), cols.size());
    A.setFromTriplets(tripletsB.begin(), tripletsB.end());
}

//
// Dense subvector.
//
// TODO: Generalize for matrices.
//
inline VecXd subDense( const VecXd& u, const std::vector<int>& rows )
{
    VecXd un = VecXd::Zero( rows.size() );
    
    for (size_t i=0; i<rows.size(); i++)
        un(i) = u(rows[i]);

    return un;
}

}   // END namespace
