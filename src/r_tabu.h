/*  Find a local maximum of a CVI via tabu-like search.
 *
 *  Copyleft (C) 2020-2021, Marek Gagolewski <https://www.gagolewski.com>
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


#ifndef R_TABU_H
#define R_TABU_H

#include <string>
#include "r_interop.h"
using namespace Rcpp;




//' (Tabu-like) (stochastic) hill climbing
//'
//' @param cvi_ptr pointer, see _CVI_create()
//' @param y0 initial label vector
//' @param allow_revisit should a `tabu` list be maintained?
//' @param max_iter_with_no_improvement how many iterations
//'        where there are no improvement in the fitness function
//'        before we give up?
//' @param max_iter maximal number of iterations
//' @param max_samples if <= 0, then an exhaustive search of all the
//'        neighbouring points is conveyed; otherwise, choose
//'        next candidates at random
//' @param verbose print additional info on the console?
//'
//' @return see optim()
//' @export
// [[Rcpp::export(".CVI_improve")]]
List _CVI_improve(
    SEXP cvi_ptr,
    NumericVector y0,
    bool allow_revisit = false,
    int max_iter_with_no_improvement = 250,
    int max_iter = 10000,
    int max_samples = -1,
    bool verbose = false)
{
    XPtr< ClusterValidityIndex > cvi =
        Rcpp::as< XPtr< ClusterValidityIndex > > (cvi_ptr);
    ClusterValidityIndex* index = &(*cvi);
    size_t K = index->get_K();
    size_t n = index->get_n();

    std::vector<uint8_t> y = translateLabels_fromR(y0);
    index->set_labels(y);

    RNGScope rngScope;


    std::unordered_set< std::vector<uint8_t>, Hash > tabuList;


    FLOAT_T best_f = index->compute();
    std::vector<uint8_t> best_y = y;
    if (!allow_revisit)
        tabuList.insert(y);

    bool random_search = true;
    if (max_samples <= 0 || max_samples >= (int)(n*K)) {
        random_search = false;
        max_samples = (int)n*K;
    }

    // bool ifChange;
    int t = 0;
    int k = 0;
    do {
        ++k;
        Rcpp::checkUserInterrupt();
        size_t  cur_best_i = 0;
        uint8_t cur_best_j = 0;
        FLOAT_T cur_best_f = -INFTY;

        // generate neighbours
        for (int s=0; s<max_samples; s++) {
            size_t i;
            uint8_t j;
            if (!random_search) {
                i = (size_t) (s/K);
                j = (uint8_t)(s%K);
            }
            else {
                i = (size_t) R::runif(0, n);
                j = (uint8_t)R::runif(0, K);
            }

            if (y[i] == j) continue;
            if (index->get_count(y[i]) <= 1) continue;

            if (!allow_revisit) {
                ssize_t tmp = y[i];
                y[i] = j;
                bool is_tabu = (tabuList.find(y) != tabuList.end());
                y[i] = tmp;
                if (is_tabu) {
                    ++t;
                    continue;
                }
            }

            index->modify(i, j);
            FLOAT_T res = index->compute();
            index->undo();

            if (res > cur_best_f){
                cur_best_f = res;
                cur_best_i = i;
                cur_best_j = j;
            }
        }

        y[cur_best_i] = cur_best_j;
        index->modify(cur_best_i, cur_best_j);

        if (!allow_revisit)
            tabuList.insert(y);

        if (cur_best_f > best_f) {
            best_f = cur_best_f;
            best_y = y;
        }
        else {
            max_iter_with_no_improvement--;
        }

        if (verbose && (k % 100 == 1 || max_iter_with_no_improvement == 0 ||
                        k == max_iter)) {
            Rprintf("%8d: best_f=%10.3f cur_best_f=%10.3f, %6d left, %6d tabus\n",
                    k, best_f, cur_best_f, max_iter_with_no_improvement, t);
        }
    }
    while (max_iter_with_no_improvement>0 && k < max_iter);

    NumericVector par = translateLabels_toR(best_y);
    return Rcpp::List::create(
        _["par"] = par,
        _["value"] = (double)best_f,
        _["counts"] = (int)k,
        _["convergence"] = 0,
        _["message"] = "max_iter_with_no_improvement or max_iter reached"
    );
}






//' Tabu-like hill climbing from multiple initial points
//'
//' an exhaustive search of all the neighbouring points is conveyed
//'
//' @param cvi_ptr pointer, see _CVI_create()
//' @param Y0 set of initial label vectors
//' @param max_iter_with_no_improvement how many iterations
//'        where there are no improvement in the fitness function
//'        before we give up?
//' @param max_iter maximal number of iterations
//' @param verbose print additional info on the console?
//'
//' @return see optim()
//' @export
// [[Rcpp::export(".CVI_improve_turbo")]]
List _CVI_improve_turbo(
    SEXP cvi_ptr,
    NumericMatrix Y0,
    int max_iter_with_no_improvement = 250,
    int max_iter = 10000,
    bool verbose = false)
{
    XPtr< ClusterValidityIndex > cvi =
        Rcpp::as< XPtr< ClusterValidityIndex > > (cvi_ptr);
    ClusterValidityIndex* index = &(*cvi);
    size_t K = index->get_K();
    size_t n = index->get_n();
    size_t max_samples = (int)n*K;
    std::unordered_set< std::vector<uint8_t>, Hash > tabuList;
    FLOAT_T best_f = -INFTY;
    std::vector<uint8_t> best_y;

    int t = 0; // tabu hits
    for (int c=0; c<Y0.ncol(); ++c) {
        std::vector<uint8_t> y = translateLabels_fromR(Y0.column(c));

        bool is_tabu = (tabuList.find(y) != tabuList.end());
        if (is_tabu) {
            ++t;
            continue;
        }


        index->set_labels(y);
        FLOAT_T cur_f = index->compute();
        if (cur_f > best_f) {
            best_f = cur_f;
            best_y = y;
        }
        tabuList.insert(y);

        if (verbose) {
            Rprintf("(%3d/%3d) %8d: best_f=%10.3f cur_best_f=%10.3f, %6d left, %6d tabu hits, %6d tabu size\r",
                c+1, Y0.ncol(), 0, best_f, cur_f, max_iter_with_no_improvement, t, tabuList.size());
        }


        // bool ifChange;
        int k = 0; // number of iterations
        int p = 0; // number of iterations with no improvement
        do {
            ++k;
            Rcpp::checkUserInterrupt();
            size_t  cur_best_i = 0;
            uint8_t cur_best_j = 0;
            FLOAT_T cur_best_f = -INFTY;

            // generate neighbours
            for (size_t s=0; s<max_samples; s++) {
                size_t i;
                uint8_t j;
                i = (size_t) (s/K);
                j = (uint8_t)(s%K);

                if (y[i] == j) continue;
                if (index->get_count(y[i]) <= 1) continue;

                ssize_t tmp = y[i];
                y[i] = j;
                bool is_tabu = (tabuList.find(y) != tabuList.end());
                y[i] = tmp;
                if (is_tabu) {
                    ++t;
                    continue;
                }

                index->modify(i, j);
                FLOAT_T res = index->compute();
                index->undo();

                if (res > cur_best_f) {
                    cur_best_f = res;
                    cur_best_i = i;
                    cur_best_j = j;
                }
            }

            if (IS_MINUS_INFTY(cur_best_f)) {
                // can't improve at all
                break;
            }

            y[cur_best_i] = cur_best_j;
            index->modify(cur_best_i, cur_best_j);

            tabuList.insert(y);

            if (cur_best_f > best_f) {
                best_f = cur_best_f;
                best_y = y;
            }
            else {
                p++;
            }

            if (verbose && (k % 10 == 1)) {
                Rprintf("(%3d/%3d) %8d: best_f=%10.3f cur_best_f=%10.3f, %6d left, %6d tabu hits, %6d tabu size\r",
                    c+1, Y0.ncol(), k, best_f, cur_best_f, max_iter_with_no_improvement-p, t, tabuList.size());
            }

        }
        while (p<max_iter_with_no_improvement && k < max_iter && !IS_PLUS_INFTY(best_f));

        if (verbose) {
            Rprintf("(%3d/%3d) %8d: best_f=%10.3f cur_best_f=%10.3f, %6d left, %6d tabu hits, %6d tabu size\r",
                c+1, Y0.ncol(), k, best_f, -INFTY, max_iter_with_no_improvement-p, t, tabuList.size());
        }

        if (IS_PLUS_INFTY(best_f)) {
            // can't improve even further
            break;
        }
    }
    if (verbose) Rprintf("\n");

    CVI_ASSERT(!IS_MINUS_INFTY(best_f)); // couldn't be worse

    NumericVector par = translateLabels_toR(best_y);
    return Rcpp::List::create(
        _["par"] = par,
        _["value"] = (double)best_f,
        //_["counts"] = (int)k,
        _["convergence"] = 0,
        _["message"] = "max_iter_with_no_improvement or max_iter reached"
    );
}


#endif
