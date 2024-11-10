#include "solvers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <optional>
#include <vector>

/**
 * @brief Gaussian elimination to solve a system of linear equations
 * @param R The matrix R in the system of linear equations
 * @param b The RHS of the system of linear equations
 * @return The solution to the system of linear equations
 * @details The function solves the system of linear equations R * x = b
 */
bool gaussian_elimination( std::vector<std::vector<int>>& R, std::vector<double>& b )
{
  uint32_t m = R.size();    /* number of rows */
  uint32_t n = R[0].size(); /* number of columns */

  for ( uint32_t i = 0; i < m; i++ )
  {
    /* Find pivot for column i */
    uint32_t max_row = i;
    for ( uint32_t j = i + 1; j < m; j++ )
      if ( std::abs( R[j][i] ) > std::abs( R[max_row][i] ) )
        max_row = j;

    /* Swap row i with max_row */
    std::swap( R[i], R[max_row] );
    std::swap( b[i], b[max_row] );

    /* Make the diagonal element 1 */
    double diag = R[i][i];
    if ( diag == 0 )
      return false;

    for ( uint32_t j = i; j < n; j++ )
      R[i][j] /= diag;
    b[i] /= diag;

    /* Make the other elements in the column 0 */
    for ( uint32_t j = 0; j < m; j++ )
    {
      if ( j == i )
        continue;
      double factor = R[j][i];
      for ( uint32_t k = i; k < n; k++ )
        R[j][k] -= factor * R[i][k];
      b[j] -= factor * b[i];
    }
  }
  for ( uint32_t i = n; i < m; i++ )
    if ( b[i] != 0 )
      return false;
  b.resize( n );
  return true;
}