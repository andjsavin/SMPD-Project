#ifndef MATRIXUTIL_H
#define MATRIXUTIL_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace bnu = boost::numeric::ublas;

int determinant_sign(const bnu::permutation_matrix<std::size_t>& pm);

double determinant( bnu::matrix<double>& m );

bnu::matrix<double> getMatrix(std::vector<std::vector<double>> v);

#endif // MATRIXUTIL_H
