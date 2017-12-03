/*
The following code inverts the matrix input using LU-decomposition with backsubstitution of unit vectors. Reference: Numerical Recipies in C, 2nd ed., by Press, Teukolsky, Vetterling & Flannery.
you can solve Ax=b using three lines of ublas code:
permutation_matrix<> piv;
lu_factorize(A, piv);
lu_substitute(A, piv, x);
*/
 #ifndef INVERT_MATRIX_HPP
 #define INVERT_MATRIX_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;
using namespace std;

 /* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
bool InvertMatrix(const matrix<double>& input, matrix<double>& inverse)
{
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<double> A(input);

    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A, pm);

    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<double> (A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}

 #endif //INVERT_MATRIX_HPP
