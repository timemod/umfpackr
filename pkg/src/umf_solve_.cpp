#include <Rcpp.h>
#include "umfpack.h"
using namespace Rcpp;

//#define TIMER
#ifdef TIMER
#include <time.h>
#endif

/*
 * This function solves x in the linear equation a  x = b
 * INPUT:
 * a        a sparse Matrix (a dgCMatrix object)
 * b        the right hand side of the equation.
 * cndtol   condition tolerance: if the rough estimate of the condition 
 *          is smaller than cndtol, then we wil not solve.
 */

// [[Rcpp::export]]
List umf_solve_(S4 a, NumericVector b, const double cndtol) {

    IntegerVector dims = a.slot("Dim");
    IntegerVector Ap = a.slot("p");
    IntegerVector Ai = a.slot("i");
    NumericVector Ax = a.slot("x");

    int n = dims[0];  // TODO: check that dims[1] == dims[0]

    NumericVector x(n);

    double *null = (double *) NULL;
    void *Symbolic, *Numeric ;
    double info[UMFPACK_INFO];

//  LU factorisation

#ifdef TIMER
    clock_t begin, end;
    begin = clock();
#endif
    (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                               &Symbolic, null, null) ;
    (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax), Symbolic, 
                               &Numeric, null, info) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    double cond = info[UMFPACK_RCOND];
#ifdef TIMER
    end = clock();
    Rprintf("Elapsed time for LU %g\n", double(end - begin) / CLOCKS_PER_SEC);
#endif

// solving
//
#ifdef TIMER
    begin = clock();
#endif
    if (cond >= cndtol) {
        (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax), 
                                 REAL(x), REAL(b), Numeric, null, null) ;
    }
    umfpack_di_free_numeric (&Numeric) ;
#ifdef TIMER
    end = clock();
    Rprintf("Elapsed time for solving %g\n", 
            double(end - begin) / CLOCKS_PER_SEC);
#endif

    if (cond >= cndtol) {
        return List::create(Named("x")    = x, Named("cond") = cond);
    } else {
        return List::create(Named("cond") = cond);
    }
}
