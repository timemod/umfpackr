#include <Rcpp.h>
#include "umfpack.h"
using namespace Rcpp;

// [[Rcpp::export]]
List umf_solve_nl_(NumericVector start, Function fn, Function jac,
                   List control) {

    bool trace = control["trace"];
    double ftol = control["ftol"];
    int maxiter = control["maxiter"];

    bool solved = false;

    int n = start.size();
    NumericVector x(n), dx(n);

    x = start;

    if (trace) {
        Rprintf("\nIteration report\n");
        Rprintf("----------------\n");
        Rprintf("%5s%20s%20s%20s\n", "Iter", "inv. cond. jac.",
                    "Largest |f|", "Index largest |f|");
    }

    double cond;
    int iter;
    NumericVector f;
    for (iter = 0; iter < maxiter; iter++) {

        f = fn(start);
        NumericVector f_abs = abs(f);
        double  err = max(f_abs);
        // TODO: handle NA values
        int i = which_max(f_abs);
        if (iter == 0) {
            Rprintf("%5d%20s%20.3e%20d\n", iter, "", err, i);
        } else {
            Rprintf("%5d%20.2e%20.3e%20d\n", iter, cond, err, i);
        }
        if (err < ftol) {
            solved = true;
            break;
        }
        S4 a = jac(start);
        IntegerVector dims = a.slot("Dim");
        IntegerVector Ap = a.slot("p");
        IntegerVector Ai = a.slot("i");
        NumericVector Ax = a.slot("x");

        // TODO: check that dims[0] == dims[1] = n?

        double *null = (double *) NULL;
        void *Symbolic, *Numeric ;

        // LU factorization
        (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                               &Symbolic, null, null) ;
        (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                   Symbolic, &Numeric, null, null) ;
        umfpack_di_free_symbolic (&Symbolic) ;

        // TODO: compute the inverse condition of the matrix
        cond = NA_REAL;

        // solve
        (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                 REAL(dx), REAL(f), Numeric, null, null) ;
        umfpack_di_free_numeric (&Numeric) ;

        x = x - dx;
    }

    return List::create(Named("solved") = solved,
                        Named("iter")   = iter,
                        Named("x")      = x,
                        Named("fval")   = f);
}
