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

#ifndef __CVI_DUNN_H
#define __CVI_DUNN_H

#include "cvi.h"




/** Dunn's index for measuring the degree to which clusters are
 *  compact and well-separated
 *
 *  The index is defined by Eq.(3) in (Dunn, 1974).
 *
 *  TODO: formula
 *
 *  TODO: time, memory complexity as a function of n, K
 *
 *  J.C. Dunn, A fuzzy relative of the ISODATA process and its use in detecting
 *  Compact Well-Separated Clusters, Journal of Cybernetics 3(3), 1974,
 *  pp. 32-57, doi:10.1080/01969727308546046.
 */
class DunnIndex : public ClusterValidityIndex
{
protected:
    matrix<DistTriple> dist; /**< intra-cluster distances:
        dist(i,j) = min( d(X(u,), X(v,)) ), X(u,) in C_i, X(v,) in C_j  (i!=j)
        */
    std::vector<DistTriple> diam; /**< cluster diameters:
        diam[i] = max( d(X(u,), X(v,)) ), X(u,), X(v,) in C_i
        */
    EuclideanDistance D; ///< squared Euclidean


    matrix<DistTriple> last_dist;      ///< for undo()
    std::vector<DistTriple> last_diam; ///< for undo()
    bool last_chg; ///< for undo() (were dist&diam changed at all?)


    void recompute_dist_diam()
    {
        for (size_t i=0; i<K; ++i) {
            diam[i] = DistTriple(0, 0, 0.0);
            //dist(i, i) = 0.0;
            for (size_t j=i+1; j<K; ++j) {
                dist(i,j) = dist(j,i) = DistTriple(0, 0, INFTY);
            }
        }

        for (size_t i=0; i<n-1; ++i) {
            for (size_t j=i+1; j<n; ++j) {
                double d = D(i, j);
                if (L[i] == L[j]) {
                    if (d > diam[L[i]].d)
                        diam[L[i]] = DistTriple(i, j, d);
                }
                else {
                    if (d < dist(L[i], L[j]).d)
                        dist(L[i], L[j]) = dist(L[j], L[i]) = DistTriple(i, j, d);
                }
            }
        }
    }


public:
    // Described in the base class
    DunnIndex(
           const matrix<FLOAT_T>& _X,
           const uint8_t _K,
           const bool _allow_undo=false)
        : ClusterValidityIndex(_X, _K, _allow_undo),
          dist(K, K),
          diam(K),
          D(&X, n<=CVI_MAX_N_PRECOMPUTE_DISTANCE, true/*squared*/),
          last_dist(K, K),
          last_diam(K)
    {

    }




    // Described in the base class
    virtual void set_labels(const std::vector<uint8_t>& _L)
    {
        ClusterValidityIndex::set_labels(_L); // sets L, count and centroids

        recompute_dist_diam();
    }


    // Described in the base class
    virtual void modify(size_t i, uint8_t j)
    {
//         uint8_t tmp = j;

        bool needs_recompute = false;
        for (size_t u=0; u<K; ++u) {
            last_diam[u] = diam[u];

            // if the point being modified determines its cluster's diameter:
            if (diam[u].i1 == i || diam[u].i2 == i)
                needs_recompute = true;

            for (size_t v=u+1; v<K; ++v) {
                // if the point being modified determines intra-cluster distance:
                if (dist(u,v).i1 == i || dist(u,v).i2 == i)
                    needs_recompute = true;

                last_dist(u,v) = last_dist(v,u) = dist(u,v);
            }
        }

        // sets L[i]=j and updates count
        ClusterValidityIndex::modify(i, j);

        if (needs_recompute) {
            last_chg = true;
            recompute_dist_diam();
        }
        else {
            last_chg = false;
            for (size_t u=0; u<n; ++u) {
                if (i == u) continue;

                double d = D(i, u);
                if (L[i] == L[u]) {
                    if (d > diam[L[i]].d) {
                        diam[L[i]] = DistTriple(i, u, d);
                        last_chg = true;
                    }
                }
                else {
                    if (d < dist(L[i], L[u]).d) {
                        dist(L[i], L[u]) = dist(L[u], L[i]) = DistTriple(i, u, d);
                        last_chg = true;
                    }
                }
            }
        }
    }


    // Described in the base class
    virtual void undo()
    {
        if (last_chg) {
            for (size_t i=0; i<K; ++i) {
                diam[i] = last_diam[i];
                for (size_t j=i+1; j<K; ++j) {
                    dist(i,j) = dist(j,i) = last_dist(i,j);
                }
            }
        }

        ClusterValidityIndex::undo();
    }


    // Described in the base class
    virtual FLOAT_T compute()
    {
        FLOAT_T max_diam = 0.0;
        FLOAT_T min_dist = INFTY;
        for (size_t i=0; i<K; ++i) {
            if (diam[i].d > max_diam)
                max_diam = diam[i].d;
            for (size_t j=i+1; j<K; ++j) {
                if (dist(i, j).d < min_dist)
                    min_dist = dist(i, j).d;
            }
        }

        return sqrt(min_dist/max_diam);
    }
};



#endif
