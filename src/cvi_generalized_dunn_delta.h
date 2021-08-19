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

#ifndef __CVI_GENERALIZED_DUNN_DELTA_H
#define __CVI_GENERALIZED_DUNN_DELTA_H

#include "cvi.h"

class Delta
{
protected:
    EuclideanDistance& D; ///< squared Euclidean
    const matrix<FLOAT_T>& X;
    //matrix<FLOAT_T>& X;         ///< data matrix of size n*d
    std::vector<uint8_t>& L;    ///< current label vector of size n
    std::vector<size_t>& count; ///< size of each of the K clusters
    uint8_t K;
    size_t n;
    size_t d;
    matrix<FLOAT_T>* centroids; ///< centroids, can be NULL
    
public:
    Delta(
           EuclideanDistance& D,
           const matrix<FLOAT_T>& X,
           std::vector<uint8_t>& L,
           std::vector<size_t>& count,
           uint8_t K,
           size_t n,
           size_t d,
           matrix<FLOAT_T>* centroids=nullptr
           )
        : D(D),
          X(X),
          L(L),
          count(count),
          K(K),
          n(n),
          d(d),
          centroids(centroids)
    { }

    virtual void before_modify(size_t i, uint8_t j) = 0;
    virtual void after_modify(size_t i, uint8_t j) = 0;
    virtual void undo() = 0;
    virtual void recompute_all() = 0;
};

class LowercaseDelta : public Delta
{
public:
    LowercaseDelta(
        EuclideanDistance& D,
        const matrix<FLOAT_T>& X,
        std::vector<uint8_t>& L,
        std::vector<size_t>& count,
        uint8_t K,
        size_t n,
        size_t d,
        matrix<FLOAT_T>* centroids=nullptr
        )
    : Delta(D,X,L,count,K,n,d,centroids)
    { }
    virtual FLOAT_T compute(size_t k, size_t l) = 0;
};


class UppercaseDelta : public Delta
{
public:
    UppercaseDelta(
        EuclideanDistance& D,
        const matrix<FLOAT_T>& X,
        std::vector<uint8_t>& L,
        std::vector<size_t>& count,
        uint8_t K,
        size_t n,
        size_t d,
        matrix<FLOAT_T>* centroids=nullptr
        )
    : Delta(D,X,L,count,K,n,d,centroids)
    { }
    virtual FLOAT_T compute(size_t k) = 0;
};

class DeltaFactory
{
public:
    virtual bool IsCentroidNeeded() = 0;
};


class LowercaseDeltaFactory : public DeltaFactory
{
public:
    // cannot be in DeltaFactory since result type is different, even if parameter list is the same
    virtual LowercaseDelta* create(EuclideanDistance& D,
           const matrix<FLOAT_T>& X,
           std::vector<uint8_t>& L,
           std::vector<size_t>& count,
           uint8_t K,
           size_t n,
           size_t d,
           matrix<FLOAT_T>* centroids=nullptr) = 0;

    // static LowercaseDeltaFactory* GetSpecializedFactory(std::string lowercaseDeltaName);
};

class UppercaseDeltaFactory : public DeltaFactory
{
public:
    // cannot be in DeltaFactory since result type is different, even if parameter list is the same
    virtual UppercaseDelta* create(EuclideanDistance& D,
           const matrix<FLOAT_T>& X,
           std::vector<uint8_t>& L,
           std::vector<size_t>& count,
           uint8_t K,
           size_t n,
           size_t d,
           matrix<FLOAT_T>* centroids=nullptr) = 0;

    // static UppercaseDeltaFactory* GetSpecializedFactory(std::string uppercaseDeltaName);
};

#endif