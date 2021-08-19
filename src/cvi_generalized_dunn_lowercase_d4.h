/*  Internal cluster validity indices
 *
 *  Copyleft (C) 2020, Maciej Bartoszuk <http://bartoszuk.rexamine.com>
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

#ifndef __CVI_GENERALIZED_DUNN_LOWERCASE_D4_H
#define __CVI_GENERALIZED_DUNN_LOWERCASE_D4_H

#include "cvi.h"
#include "cvi_generalized_dunn_delta.h"

class LowercaseDelta4 : public LowercaseDelta
{
public:
    LowercaseDelta4(
        EuclideanDistance& D,
        const matrix<FLOAT_T>& X,
        std::vector<uint8_t>& L,
        std::vector<size_t>& count,
        uint8_t K,
        size_t n,
        size_t d,
        matrix<FLOAT_T>* centroids=nullptr
        )
    : LowercaseDelta(D, X, L, count,K,n,d,centroids)
    { 
    }
    virtual void before_modify(size_t i, uint8_t j) {
        // all happens in CentroidsBasedIndex
    }
    virtual void after_modify(size_t i, uint8_t j) {
        // all happens in CentroidsBasedIndex
    }
    virtual void undo() {
        // all happens in CentroidsBasedIndex
    }
    virtual void recompute_all() {
        // all happens in CentroidsBasedIndex
    }
    
    virtual FLOAT_T compute(size_t k, size_t l) {
        FLOAT_T act = 0.0;
        for (size_t u=0; u<d; ++u) {
            act += square((*centroids)(k, u) - (*centroids)(l, u));
        }
        return sqrt(act);
    }

}; 

class LowercaseDelta4Factory : public LowercaseDeltaFactory
{
public:
    virtual bool IsCentroidNeeded() { return true; }
    
    virtual LowercaseDelta* create(EuclideanDistance& D,
           const matrix<FLOAT_T>& X,
           std::vector<uint8_t>& L,
           std::vector<size_t>& count,
           uint8_t K,
           size_t n,
           size_t d,
           matrix<FLOAT_T>* centroids=nullptr) {
               return new LowercaseDelta4(D, X, L, count, K, n, d, centroids);
           }
};

#endif