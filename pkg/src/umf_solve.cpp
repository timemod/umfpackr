#include <Rcpp.h>
#include "umfpack.h"
using namespace Rcpp;

//#define TIMER
#ifdef TIMER
#include <time.h>
#endif

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

//  LU factorisation
#ifdef TIMER
    clock_t begin, end;
    begin = clock();
#endif
    (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                               &Symbolic, null, null) ;
    (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax), Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic) ;
#ifdef TIMER
    end = clock();
    Rprintf("Elapsed time for LU %g\n", double(end - begin) / CLOCKS_PER_SEC);
#endif

// solving
#ifdef TIMER
    begin = clock();
#endif
    (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax), 
            REAL(x), REAL(b), Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;
#ifdef TIMER
    end = clock();
    Rprintf("Elapsed time for solving %g\n", 
            double(end - begin) / CLOCKS_PER_SEC);
#endif

    return x;
}
