/*  Internal cluster validity indices
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

#ifndef __CVI_GAMMA_H
#define __CVI_GAMMA_H

#include "cvi.h"





/** The Baker-Hubert Gamma Coefficient
 *
 *  Equal to the Goodman-Kruskal rank correlation measure
 *  between the vector of all pairwise distances and the 0/1 indicator
 *  function stating whether the corresponding pair of points belongs
 *  to the same cluster (0) or not (1).
 *
 *  Gives a value between -1 and 1. Gamma=(NC-ND)/(NC+ND),
 *  where NC - number of concordant and ND - number of discordant pairs.
 *
 *  TODO: formula
 *
 *  Time complexity: O(n^2).
 *  Memory complexity: O(n^2).
 *
 *  F.B. Baker, L.J. Hubert, Measuring the power of hierarchical cluster
 *  analysis, Journal of the American Statistical Association 70(349), 1975,
 *  pp. 31-38.
 *
 */
class GammaIndex : public ClusterValidityIndex
{
protected:
    size_t n_pairs; ///< n*(n-1)/2
    std::vector<DistTriple> D;

public:
    // Described in the base class
    GammaIndex(
           const matrix<FLOAT_T>& _X,
           const uint8_t _K,
           const bool _allow_undo=false)
        : ClusterValidityIndex(_X, _K, _allow_undo),
            n_pairs(n*(n-1)/2),
            D(n_pairs)

    {
        size_t k=0;
        for (size_t i=0; i<n-1; ++i) {
            for (size_t j=i+1; j<n; ++j) {
                D[k++] = DistTriple(i, j,
                    distance_l2_squared(X.row(i), X.row(j), X.ncol()));
            }
        }
        std::sort(D.begin(), D.end());
    }


    // Described in the base class
    virtual FLOAT_T compute()
    {
        // let I1, I2 - an auxiliary vector of size n_pairs=n*(n-1)/2
        // with (I1[i], I2[i]) - indexes of the i-th nearest pair of points,
        // i.e., D(I1[0], I2[0]) <= ... D(I1[i], I2[i]) ... <= D(I1[n_pairs], I2[n_pairs])


        size_t nc = 0; ///< number of concordant pairs
        size_t nd = 0; ///< number of discordant pairs

        size_t number_of_0_so_far = 0;
        size_t number_of_1_so_far = 0;
        for (size_t j=0; j<n_pairs; ++j) {
            if (L[D[j].i1] == L[D[j].i2]) { // 0 - a pair of objects in the same cluster
                nd += number_of_1_so_far;
                number_of_0_so_far++;
            }
            else { // 1
                nc += number_of_0_so_far;
                number_of_1_so_far++;
            }
        }

        FLOAT_T ret = ((FLOAT_T)nc-(FLOAT_T)nd)/((FLOAT_T)nc+(FLOAT_T)nd);
        CVI_ASSERT(std::fabs(ret) < 1.0+1e-9);
        return ret;
    }
};



#endif
