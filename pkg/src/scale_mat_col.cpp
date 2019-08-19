#include <Rcpp.h>
#include <cmath>
#include "umfpack.h"
#include <algorithm>
using namespace Rcpp;

// This function applies column scaling to matrix mat.
// and returns the scaling factors.
// Note that this function has a side effect: the input
// arguments are modified.
// INPUT AND OUTPUT (the input values are modified):
//    mat    a matrix in compressed column format
//    scale  a list with dimension n (size of the matrix).
//  RETURNS
//   scale
// [[Rcpp::export]]
NumericVector scale_mat_col(S4 mat, NumericVector scale) {

    // This implementation is similar to the implementation in 
    // nlqslv.

    IntegerVector dims = mat.slot("Dim");
    IntegerVector Ap = mat.slot("p");
    IntegerVector Ai = mat.slot("i");
    NumericVector Ax = mat.slot("x");

    int n = dims[0];

    for (int c = 0; c < n; c++) {
        int i1 = Ap[c];
        int i2 = Ap[c + 1];
        double s = 0;
        for (int i = i1; i < i2; i++) {
            double value = Ax[i];
            if (ISNAN(value)) {
                s = value;
                break;
            }
            s = s + std::abs(value);
        }
        if (ISNAN(s) || s == 0) {
            s = 1;
        }
        scale[c] = s;
    }

    // scale the matrix
    for (int c = 0; c < n; c++) {
        int i1 = Ap[c];
        int i2 = Ap[c + 1];
        for (int i = i1; i < i2; i++) {
            Ax[i] = Ax[i] / scale[c];
        }
    }

    return scale;
}
