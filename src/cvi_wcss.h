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

#ifndef __CVI_WCSS_H
#define __CVI_WCSS_H

#include "cvi.h"



/** Negated Within-Cluster Sum of Squares and the Ball-Hall Index
 *
 *  The Ball-Hall index is weighted by the cluster cardinality (weighted=true).
 *
 *  TODO: formula
 *
 *  TODO: time, memory complexity as a function of n, K
 *
 *  WCSS is the objective function used, amongst others, in the k-means and
 *  the Ward and Caliński&Harabasz algorithms.
 *  The Ball-Hall index
 *
 *  TODO: cite Ward;Edwards&Cavalli-Sforza;
 *  MacQueen, 1967; LLoyd, 1957
 *
 *  G.H. Ball, D.J. Hall,
 *  ISODATA: A novel method of data analysis and pattern classification,
 *  Technical report No. AD699616, Stanford Research Institute, 1965.
 *
 *  T. Caliński, J. Harabasz, A dendrite method for cluster analysis,
 *  Communications in Statistics, 3(1), 1974, pp. 1-27,
 *  doi:10.1080/03610927408827101.
 *
 */
class WCSSIndex : public CentroidsBasedIndex
{
protected:
    bool weighted;          ///< false for WCSS, true for the Ball-Hall index

public:
    // Described in the base class
    WCSSIndex(
           const matrix<FLOAT_T>& _X,
           const uint8_t _K,
           const bool _allow_undo=false,
           bool _weighted=false)
        : CentroidsBasedIndex(_X, _K, _allow_undo)
    {
        weighted = _weighted;
    }



    // Described in the base class
    virtual FLOAT_T compute()
    {
        // sum of within-cluster squared L2 distances
        FLOAT_T wcss = 0.0;
        for (size_t i=0; i<n; ++i) {
            for (size_t j=0; j<d; ++j) {
                wcss += square(centroids(L[i],j)-X(i,j))/((weighted)?count[L[i]]:1.0);
            }
        }
        return -wcss;  // negative!!!
    }

};



#endif
