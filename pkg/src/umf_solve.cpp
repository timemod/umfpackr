#include <Rcpp.h>
#include "umfpack.h"
using namespace Rcpp;


//' Solve linear equation using UMFPACK
//'
//' @param A a \code{dgCMatrix}
//' @param b the right hand side of A x = b
//' @return x the solution x.
//' @export
//' @useDynLib umfpackr
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericVector umf_solve(S4 a, NumericVector b) {

    IntegerVector dims = a.slot("Dim");
    IntegerVector Ap = a.slot("p");
    IntegerVector Ai = a.slot("i");
    NumericVector Ax = a.slot("x");

    int n = dims[0];  // TODO: check that dims[1] == dims[0]

    NumericVector x(n);

    double *null = (double *) NULL;
    void *Symbolic, *Numeric ;
    (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                               &Symbolic, null, null) ;
    (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax), Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax), 
            REAL(x), REAL(b), Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;

    return x;
}
