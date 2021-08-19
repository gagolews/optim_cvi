/*  Common functions, macros, includes
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


#ifndef __c_common_h
#define __c_common_h

#include <stdexcept>
#include <string>
#include <limits>
#include <vector>

typedef double FLOAT_T; ///< float type we are working internally with


template <typename T> inline T square(const T x) { return x*x; }





/** Class to hash a vector of small ints.
 */
struct Hash {
    /** Computes a hash of a given vector of small ints.
     *
     * Based on the boost library's hash_combine()
     *
     * @param x
     * @return
     */
    size_t operator() (const std::vector<uint8_t>& x) const {
        size_t seed = x[0];
        size_t n = x.size();
        for (size_t i=0; i<n; ++i) {
            // boost::hash_combine()
            seed ^= (size_t)x[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};



#define CVI_MAX_N_PRECOMPUTE_DISTANCE 10000

#ifndef CVI_ASSERT
#define __CVI_STR(x) #x
#define CVI_STR(x) __CVI_STR(x)

#define CVI_ASSERT(EXPR) { if (!(EXPR)) \
    throw std::runtime_error( "CVI: Assertion " #EXPR " failed in "\
        __FILE__ ":" CVI_STR(__LINE__) ); }
#endif

#ifndef GENIECLUST_ASSERT
#define GENIECLUST_ASSERT CVI_ASSERT
#endif


#ifndef INFTY
#define INFTY (std::numeric_limits<FLOAT_T>::infinity())
#endif

#define IS_PLUS_INFTY(x)  ((x) > 0.0 && !std::isfinite(x))
#define IS_MINUS_INFTY(x) ((x) < 0.0 && !std::isfinite(x))

#endif
