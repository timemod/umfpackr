#include <Rcpp.h>
#include <umfpack.h>
#include "set_umf_control.h"
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
 */

// [[Rcpp::export]]
List umf_solve_(S4 a, NumericVector b, List umf_control) {

    IntegerVector dims = a.slot("Dim");
    IntegerVector Ap = a.slot("p");
    IntegerVector Ai = a.slot("i");
    NumericVector Ax = a.slot("x");

    int n = dims[0];  // TODO: check that dims[1] == dims[0]

    NumericVector x(n);

    void *Symbolic, *Numeric ;
    double info[UMFPACK_INFO], control[UMFPACK_CONTROL];

    set_umf_control(control, umf_control);

//
//  LU factorisation
//  

#ifdef TIMER
    clock_t begin, end;
    begin = clock();
#endif

    (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                               &Symbolic, control, info) ;
    if (info[UMFPACK_STATUS] != UMFPACK_OK) {
        umfpack_di_free_symbolic(&Symbolic);
        if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory) {
            Rf_error("Not enough memory for symbolic factorization");
        } else {
            Rf_error("Unknown error in symbolic factorization");
        }
        return R_NilValue;
    }

    (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax), Symbolic, 
                               &Numeric, control, info) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    double stat = info[UMFPACK_STATUS]; 
    double cond = info[UMFPACK_RCOND];
    if (stat == UMFPACK_WARNING_singular_matrix) {
      return List::create(Named("status") = "singular matrix", 
                          Named("cond") = info[UMFPACK_RCOND]);
    } else if (stat != UMFPACK_OK) {
        umfpack_di_free_numeric(&Numeric);
        if (stat == UMFPACK_ERROR_out_of_memory) {
            Rf_error("Not enough memory for numeric factorization");
        } else {
            Rf_error("Unknown error in numeric factorization");
        }
        return R_NilValue;
    }
#ifdef TIMER
    end = clock();
    Rprintf("Elapsed time for LU %g\n", double(end - begin) / CLOCKS_PER_SEC);
#endif

//
// solving
//

#ifdef TIMER
    begin = clock();
#endif
   (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax), 
                             REAL(x), REAL(b), Numeric, control, info) ;
    if (info[UMFPACK_STATUS] != UMFPACK_OK) {
        umfpack_di_free_numeric (&Numeric) ;
        if (info[UMFPACK_STATUS] == UMFPACK_ERROR_out_of_memory) {
            Rf_error("Not enough memory to solve to linear system");
        } else {
            Rf_error("Unknown error while solving the linear system");
        }
        return R_NilValue;
    }

    umfpack_di_free_numeric (&Numeric) ;

#ifdef TIMER
    end = clock();
    Rprintf("Elapsed time for solving %g\n", 
            double(end - begin) / CLOCKS_PER_SEC);
#endif
    return List::create(Named("status") = "OK", Named("x") = x, 
                        Named("cond") = cond);
}
