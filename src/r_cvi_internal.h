/*  Rcpp exports - CVI
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


#ifndef R_CVI_INTERNAL_H
#define R_CVI_INTERNAL_H

#include <string>
#include <cstring>
#include "r_interop.h"
#include "cvi_calinski_harabasz.h"
#include "cvi_davies_bouldin.h"
#include "cvi_silhouette.h"
#include "cvi_dunn.h"
#include "cvi_wcss.h"
#include "cvi_wcnn.h"
#include "cvi_dunnowa.h"
#include "cvi_gamma.h"
#include "cvi_generalized_dunn.h"
#include "cvi_generalized_dunn_lowercase_d1.h"
#include "cvi_generalized_dunn_lowercase_d2.h"
#include "cvi_generalized_dunn_lowercase_d3.h"
#include "cvi_generalized_dunn_lowercase_d4.h"
#include "cvi_generalized_dunn_lowercase_d5.h"
#include "cvi_generalized_dunn_lowercase_d6.h"
#include "cvi_generalized_dunn_uppercase_d1.h"
#include "cvi_generalized_dunn_uppercase_d2.h"
#include "cvi_generalized_dunn_uppercase_d3.h"
// #include "cvi_davies_bouldin_star.h"
// #include "cvi_score_function.h"
// #include "cvi_cs.h"
// #include "cvi_cop.h"
// #include "cvi_sym.h"
using namespace Rcpp;



//' @export
// [[Rcpp::export(".CVI_create")]]
SEXP _CVI_create(Rcpp::String type, NumericMatrix X, int K, bool allow_undo=true)
{
    ClusterValidityIndex* cvi;

    const char* _type = type.get_cstring();

    if (type == "CalinskiHarabasz") {
        cvi = new CalinskiHarabaszIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo);
    }
    else if (type == "DaviesBouldin") {
        cvi = new DaviesBouldinIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo);
    }
    else if (type == "Silhouette") {
        cvi = new SilhouetteIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo, false);
    }
    else if (type == "SilhouetteW") {
        cvi = new SilhouetteIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo, true);
    }
    else if (type == "Dunn") {
        cvi = new DunnIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo);
    }
    else if (type == "WCSS") {
        cvi = new WCSSIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo, false/*not weighted*/);
    }
    else if (type == "BallHall") {
        cvi = new WCSSIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo, true/*weighted*/);
    }
    else if (type == "Gamma") {
        cvi = new GammaIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo);
    }
    else if (strncmp(_type, "DuNN_", 5) == 0) { // DuNN_M_numerator_denominator
        // e.g., DuNN_25_Min_Max
        int M = 10;
        int owa_numerator = OWA_MIN;
        int owa_denominator = OWA_MAX;

        std::string opts = std::string(_type+5);
        size_t pos1 = opts.find('_');
        CVI_ASSERT(pos1 != std::string::npos);
        M = std::atoi(opts.substr(0, pos1).c_str());
        CVI_ASSERT(M>0);  // M = min(n-1, M) in the constructor

        size_t pos2 = opts.find('_', pos1+1);
        CVI_ASSERT(pos2 != std::string::npos);
        std::string owa_numerator_str = opts.substr(pos1+1, pos2-pos1-1);
        std::string owa_denominator_str = opts.substr(pos2+1);

        owa_numerator = DuNNOWA_get_OWA(owa_numerator_str);
        owa_denominator = DuNNOWA_get_OWA(owa_denominator_str);

        cvi = new DuNNOWAIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo, M, owa_numerator, owa_denominator);
    }
    else if (strncmp(_type, "WCNN_", 5) == 0) { // WCNN_M
        int M = 0;
        if (_type[5] >= '0' && _type[5] <= '9')
            M = std::atoi(_type+5);
        CVI_ASSERT(M>0);  // M = min(n-1, M) in the constructor

        cvi = new WCNNIndex(
            matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
            K, allow_undo, M);
    }
    else if (strncmp(_type, "GDunn_", 6) == 0) {
        std::string type_string = std::string(_type);
        std::string numeratorDeltaName = type_string.substr(6, 2);
        std::string denominatorDeltaName = type_string.substr(9, 2);
        bool doesNumeratorNeedCentroids;
        bool doesDenominatorNeedCentroids;
        LowercaseDeltaFactory* lowercaseDeltaFactory;
        UppercaseDeltaFactory* uppercaseDeltaFactory;

        //Rcpp::Rcout << numeratorDeltaName << " " << denominatorDeltaName << std::endl;
        //lowercaseDeltaFactory = LowercaseDeltaFactory::GetSpecializedFactory(numeratorDeltaName);

        if (numeratorDeltaName == "d1") {
            lowercaseDeltaFactory = new LowercaseDelta1Factory();
        }
        else if (numeratorDeltaName == "d2") {
            lowercaseDeltaFactory = new LowercaseDelta2Factory();
        }
        else if (numeratorDeltaName == "d3") {
            lowercaseDeltaFactory = new LowercaseDelta3Factory();
        }
        else if (numeratorDeltaName == "d4") {
            lowercaseDeltaFactory = new LowercaseDelta4Factory();
        }
        else if (numeratorDeltaName == "d5") {
            lowercaseDeltaFactory = new LowercaseDelta5Factory();
        }
        else if (numeratorDeltaName == "d6") {
            lowercaseDeltaFactory = new LowercaseDelta6Factory();
        }
        else {
            Rf_error("invalid numeratorDeltaName (d?)");
        }

        // uppercaseDeltaFactory = UppercaseDeltaFactory::GetSpecializedFactory(denominatorDeltaName);
        if (denominatorDeltaName == "D1") {
            uppercaseDeltaFactory = new UppercaseDelta1Factory();
        }
        else if (denominatorDeltaName == "D2") {
            uppercaseDeltaFactory = new UppercaseDelta2Factory();
        }
        else if (denominatorDeltaName == "D3") {
            uppercaseDeltaFactory = new UppercaseDelta3Factory();
        }
        else {
            Rf_error("invalid denominatorDeltaName (D?)");
        }

        bool areCentroidsNeeded = lowercaseDeltaFactory->IsCentroidNeeded() || uppercaseDeltaFactory->IsCentroidNeeded();
        if (areCentroidsNeeded) {
            cvi = new GeneralizedDunnIndexCentroidBased(
                matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
                K,
                lowercaseDeltaFactory,
                uppercaseDeltaFactory,
                allow_undo);
        }
        else {
            cvi = new GeneralizedDunnIndex(
                matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
                K,
                lowercaseDeltaFactory,
                uppercaseDeltaFactory,
                allow_undo);
        }
        delete lowercaseDeltaFactory;
        delete uppercaseDeltaFactory;
    }
    else {
        Rf_error("invalid type");
    }
//     } else if (type == "Sym") {
//         cvi = new SymIndex(
//             matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
//             K, allow_undo);
//     } else if (type == "CS") {
//         cvi = new CSIndex(
//             matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
//             K, allow_undo);
//     } else if (type == "COP") {
//         cvi = new COPIndex(
//             matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
//             K, allow_undo);
//     } else if (type == "DaviesBouldinStar") {
//         cvi = new DaviesBouldinStarIndex(
//             matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
//             K, allow_undo);
//     } else if (type == "ScoreFunction") {
//         cvi = new ScoreFunctionIndex(
//             matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
//             K, allow_undo);
//     } else if (type == "SymDB") {
//         cvi = new SymDBIndex(
//             matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
//             K, allow_undo);

    XPtr< ClusterValidityIndex > retval =
        XPtr< ClusterValidityIndex >((ClusterValidityIndex*)cvi, true);

    return retval;
}


//' @export
// [[Rcpp::export(".CVI_set_labels")]]
void _CVI_set_labels(SEXP cvi_ptr, NumericVector y)
{
    XPtr< ClusterValidityIndex > cvi =
        Rcpp::as< XPtr< ClusterValidityIndex > > (cvi_ptr);
    (*cvi).set_labels(translateLabels_fromR(y));
}


//' @export
// [[Rcpp::export(".CVI_compute")]]
double _CVI_compute(SEXP cvi_ptr)
{
    XPtr< ClusterValidityIndex > cvi =
        Rcpp::as< XPtr< ClusterValidityIndex > > (cvi_ptr);
    return (double)((*cvi).compute());
}



//' @export
// [[Rcpp::export(".CVI_undo")]]
void _CVI_undo(SEXP cvi_ptr)
{
    XPtr< ClusterValidityIndex > cvi =
        Rcpp::as< XPtr< ClusterValidityIndex > > (cvi_ptr);
    (*cvi).undo();
}


//' @export
// [[Rcpp::export(".CVI_modify")]]
void _CVI_modify(SEXP cvi_ptr, int i, int j)
{ // uses 1-based indexing
    XPtr< ClusterValidityIndex > cvi =
        Rcpp::as< XPtr< ClusterValidityIndex > > (cvi_ptr);
    (*cvi).modify(i-1, j-1);
}











//' @title The Calinski-Harabasz Cluster Validity Index (Variance Ratio Criterion)
//'
//' TODO: update this docstring
//'
//' @references
//' T. Calinski & J. Harabasz. A dendrite method for cluster analysis,
//' Communications in Statistics, 3(1), 1974, pp. 1-27,
//' doi:10.1080/03610927408827101.
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//' @export
// [[Rcpp::export]]
double CVI_CalinskiHarabasz(NumericMatrix X, NumericVector y, int K)
{
    CalinskiHarabaszIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}


//' @title Negated Within-Cluster Sum of Squares
//'
//' TODO: update this docstring
//'
//' Objective function used, amongst others, in the k-means and
//' the Ward and Caliński&Harabasz algorithms.
//'
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//' @export
// [[Rcpp::export]]
double CVI_WCSS(NumericMatrix X, NumericVector y, int K)
{
    WCSSIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K, false, false/*not weighted*/
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}


//' @title Negated Ball-Hall Index
//'
//' TODO: update this docstring
//'
//' Within cluster sum of squares weighted by the cluster cardinality
//'
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//' @export
// [[Rcpp::export]]
double CVI_BallHall(NumericMatrix X, NumericVector y, int K)
{
    WCSSIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K, false, true/*weighted*/
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}



//' @title The Baker-Hubert Gamma Coefficient
//'
//' TODO: update this docstring
//'
//' Equal to the Goodman-Kruskal rank correlation measure
//' between the vector of all pairwise distances and the 0/1 indicator
//' function stating whether the corresponding pair of points belongs
//' to the same cluster (0) or not (1).
//'
//' Gives a value between -1 and 1. (NC-ND)/(NC+ND),
//' where NC - number of concordant and ND - number of discordant pairs.
//'
//' F.B. Baker, L.J. Hubert, Measuring the power of hierarchical cluster
//' analysis, Journal of the American Statistical Association 70(349), 1975,
//' pp. 31-38.
//'
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//' @export
// [[Rcpp::export]]
double CVI_Gamma(NumericMatrix X, NumericVector y, int K)
{
    GammaIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}



//' @title The Negated Davies-Bouldin Cluster Validity Index
//'
//' TODO: update this docstring
//'
//' @references
//' D.L. Davies, D.W. Bouldin,
//' A Cluster Separation Measure,
//' IEEE Transactions on Pattern Analysis and Machine Intelligence. PAMI-1 (2),
//' 1979, pp. 224-227, doi:10.1109/TPAMI.1979.4766909
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index, which is the additive inverse (negation)
//'         of the original Davies-Bouldin score. Negation is for consistency
//'         with other indices -- the higher the score, the better.
//' @export
// [[Rcpp::export]]
double CVI_DaviesBouldin(NumericMatrix X, NumericVector y, int K)
{
    DaviesBouldinIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}





//' @title The Silhouette Coefficient
//'
//' TODO: update this docstring
//'
//' @references
//' P.J. Rousseeuw, Silhouettes: a Graphical Aid to the Interpretation and
//' Validation of Cluster Analysis, Computational and Applied Mathematics 20,
//' 1987, pp. 53-65, doi:10.1016/0377-0427(87)90125-7.
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//'
//' @export
// [[Rcpp::export]]
double CVI_Silhouette(NumericMatrix X, NumericVector y, int K)
{
    SilhouetteIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K, false, false
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}


//' @title The Silhouette Coefficient
//'
//' TODO: update this docstring
//'
//' @references
//' P.J. Rousseeuw, Silhouettes: a Graphical Aid to the Interpretation and
//' Validation of Cluster Analysis, Computational and Applied Mathematics 20,
//' 1987, pp. 53-65, doi:10.1016/0377-0427(87)90125-7.
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//'
//' @export
// [[Rcpp::export]]
double CVI_SilhouetteW(NumericMatrix X, NumericVector y, int K)
{
    SilhouetteIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K, false, true
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}


//' @title Dunn's index for measuring the degree to which clusters are compact separated
//'
//' TODO: update this docstring
//'
//' The index is defined by Eq.(3) in (Dunn, 1973).
//'
//' @references
//' J.C. Dunn, A Fuzzy Relative of the ISODATA Process and Its Use in Detecting
//' Compact Well-Separated Clusters, Journal of Cybernetics 3(3), 1973,
//' pp. 32-57, doi:10.1080/01969727308546046.
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//'
//' @export
// [[Rcpp::export]]
double CVI_Dunn(NumericMatrix X, NumericVector y, int K)
{
    DunnIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K
    );
    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}


//' @title Generalised Dunn's index
//'
//' TODO: update this docstring
//'
//' @references
//' TODO: update
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//'
//' @return The computed index.
//'
//' @export
// [[Rcpp::export]]
double CVI_GDunn(NumericMatrix X, NumericVector y, int K, int lowercaseDelta, int uppercaseDelta)
{
    LowercaseDeltaFactory* lowercaseDeltaFactory;
    UppercaseDeltaFactory* uppercaseDeltaFactory;

    if (lowercaseDelta == 1) {
        lowercaseDeltaFactory = new LowercaseDelta1Factory();
    }
    else if (lowercaseDelta == 2) {
        lowercaseDeltaFactory = new LowercaseDelta2Factory();
    }
    else if (lowercaseDelta == 3) {
        lowercaseDeltaFactory = new LowercaseDelta3Factory();
    }
    else if (lowercaseDelta == 4) {
        lowercaseDeltaFactory = new LowercaseDelta4Factory();
    }
    else if (lowercaseDelta == 5) {
        lowercaseDeltaFactory = new LowercaseDelta5Factory();
    }
    else if (lowercaseDelta == 6) {
        lowercaseDeltaFactory = new LowercaseDelta6Factory();
    }
    else {
        Rf_error("invalid lowercaseDelta");
    }

    if (uppercaseDelta == 1) {
        uppercaseDeltaFactory = new UppercaseDelta1Factory();
    }
    else if (uppercaseDelta == 2) {
        uppercaseDeltaFactory = new UppercaseDelta2Factory();
    }
    else if (uppercaseDelta == 3) {
        uppercaseDeltaFactory = new UppercaseDelta3Factory();
    }
    else {
        Rf_error("invalid uppercaseDelta");
    }

    bool areCentroidsNeeded = lowercaseDeltaFactory->IsCentroidNeeded() || uppercaseDeltaFactory->IsCentroidNeeded();
    if (areCentroidsNeeded) {
        GeneralizedDunnIndexCentroidBased ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        K,
        lowercaseDeltaFactory,
        uppercaseDeltaFactory
        );

        delete lowercaseDeltaFactory;
        delete uppercaseDeltaFactory;

        ind.set_labels(translateLabels_fromR(y));
        return (double)ind.compute();
    } else {
        GeneralizedDunnIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K,
        lowercaseDeltaFactory,
        uppercaseDeltaFactory
        );

        delete lowercaseDeltaFactory;
        delete uppercaseDeltaFactory;

        ind.set_labels(translateLabels_fromR(y));
        return (double)ind.compute();
    }
}


//' @title Within-Cluster Nearest-Neighbours
//'
//' TODO: update this docstring
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//' @param M number of nearest neighbours
//'
//' @return The computed index.
//'
//' @export
// [[Rcpp::export]]
double CVI_WCNN(NumericMatrix X, NumericVector y, int K, int M=10)
{
    CVI_ASSERT(M>0);  // M = min(n-1, M) in the constructor

    WCNNIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K, false, M
    );

    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}



//' @title OWA-based Dunn-like Indices Based on Near Neighbours
//'
//' TODO: update this docstring
//'
//' @param X data matrix of size n*d
//' @param y vector of n integer labels in [1, K], where `y[i]`
//'          is the cluster id of the i-th point, `X[i,]`
//' @param K number of clusters, `max(y)`
//' @param M number of nearest neighbours
//'
//' @return The computed index.
//'
//' @export
// [[Rcpp::export]]
double CVI_DuNNOWA(NumericMatrix X, NumericVector y, int K, int M=10,
                Rcpp::String owa_numerator="Min",
                Rcpp::String owa_denominator="Max")
{
    CVI_ASSERT(M>0);    // M = min(n-1, M) in the constructor

    int _owa_numerator = DuNNOWA_get_OWA(std::string(owa_numerator));
    int _owa_denominator = DuNNOWA_get_OWA(std::string(owa_denominator));

    DuNNOWAIndex ind(
        matrix<FLOAT_T>(REAL(SEXP(X)), X.nrow(), X.ncol(), false),
        (uint8_t)K, false, M, _owa_numerator, _owa_denominator
    );

    ind.set_labels(translateLabels_fromR(y));
    return (double)ind.compute();
}



#endif
