/*  Functions for the with-R interoperability
 *
 *  Copyleft (C) 2020, Marek Gagolewski <https://www.gagolewski.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License
 *  Version 3, 19 November 2007, published by the Free Software Foundation.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Affero General Public License Version 3 for more details.
 *  You should have received a copy of the License along with this program.
 *  If this is not the case, refer to <https://www.gnu.org/licenses/>.
 */



#ifndef __R_INTEROP_H
#define __R_INTEROP_H

#include "common.h"
#include "matrix.h"
#include "cvi.h"
#include <Rcpp.h>
#include <vector>



/** Converts a 1-based label vector to a 0-based vector of small integers.
 *
 * @param x numeric vector with integer elements in [1, 256].
 * @return
 */
std::vector<uint8_t> translateLabels_fromR(const Rcpp::NumericVector& x)
{
    size_t n = x.size();
    std::vector<uint8_t> ret(n);
    for (size_t i=0; i<n; ++i) {
        int xi = (int)x[i];
        CVI_ASSERT(xi >= 1 && xi <= 256)
        ret[i] = (uint8_t)(xi-1); // 1-based -> 0-based
    }
    return ret;
}


/** Converts a 0-based label vector to a 1-based one.
 *
 * @param x
 * @return
 */
Rcpp::NumericVector translateLabels_toR(const std::vector<uint8_t>& x)
{
    size_t n = x.size();
    Rcpp::NumericVector ret(n);
    for (size_t i=0; i<n; ++i)
        ret[i] = (x[i]+1); // 1-based -> 0-based
    return ret;
}


/** Convert Rcpp's numeric matrix object (column-major) to our
 * internal type (row-major).
 *
 * @param X
 * @return
 */
matrix<FLOAT_T> translateMatrix_fromR(const Rcpp::NumericMatrix& X)
{
//     size_t n = X.nrow();
//     size_t d = X.ncol();
//     matrix<FLOAT_T> Y(n, d);
//     for (size_t i=0; i<n; i++) {
//         for (size_t j=0; j<d; j++) {
//                 Y(i, j) = X(i, j);
//         }
//     }
//     return Y;
    return matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false);
}


#endif
