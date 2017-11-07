#include "matrixutil.hpp"

int determinant_sign(const bnu::permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}

double determinant(bnu::matrix<double> &m ) {
    bnu::permutation_matrix<std::size_t> pm(m.size1());
    double det = 1.0;
    if( bnu::lu_factorize(m,pm) ) {
        det = 0.0;
    } else {
        for(int i = 0; i < m.size1(); i++)
            det *= m(i,i); // multiply by elements on diagonal
        det = det * determinant_sign( pm );
    }
    return det;
}

bnu::matrix<double> getMatrix(std::vector<std::vector<double>> v)
{
    bnu::matrix<double> m (v.size(), v.size());
    for (int i = 0; i < m.size1(); i++)
        for (int j = 0; j < m.size2(); j++) {
            m(i, j) = v[i][j];
        }
    return m;
}
