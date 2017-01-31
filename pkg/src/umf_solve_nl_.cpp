#include <Rcpp.h>
#include "umfpack.h"
#include <time.h>

using namespace Rcpp;

clock_t t_begin;
double t_f, t_jac, t_lu, t_solve;

// [[Rcpp::export]]
List umf_solve_nl_(NumericVector start, Function fn, Function jac,
                   List control) {

    bool trace    = control["trace"];
    double ftol   = control["ftol"];
    int maxiter   = control["maxiter"];

    //initialize initial convergence speed threshold

    bool solved = false;

    if (trace) {
        Rprintf("\nIteration report\n");
        Rprintf("----------------\n");
        Rprintf("%5s%20s%20s%20s\n", "Iter", "inv. cond. jac.",
                    "Largest |f|", "Index largest |f|");
    }

    int iter;

    double *null = (double *) NULL;
    double info[UMFPACK_INFO];
    double cond = NA_REAL;

    int n = start.size();
    NumericVector f;
    NumericVector x(clone(start));
    double *dx = new double[n];

    t_f     = 0.0;
    t_jac   = 0.0;
    t_lu    = 0.0;
    t_solve = 0.0;

    for (iter = 0; iter < maxiter; iter++) {

        t_begin = clock();
        f = fn(x);
        t_f += double(clock() - t_begin) / CLOCKS_PER_SEC;

        NumericVector f_abs = abs(f);
        double fmax = max(f_abs);

        // TODO: handle NA values
        if (R_IsNA(fmax)) {
            int i_na = NA_INTEGER;
            for (int i = 0; i < n; i++) {
                if (R_IsNA(f[i])) {
                    i_na = i + 1;
                    break;
                }                
            }
            if (iter == 0) {
                Rprintf("Initial value of function contains " 
                       "non-finite values (starting at index=%d)\n", i_na);
            } else {
                Rprintf("Function value contains "
                    "non-finite values (starting at index=%d)\n", i_na);
            }
            break;
        }

        if (trace) {
            int i = which_max(f_abs);
            if (iter > 0) {
                Rprintf("%5d%20.2e%20.3e%20d\n", iter, cond, fmax, i);
            } else {
                Rprintf("%5d%20s%20.3e%20d\n", iter, "", fmax, i);
            }
        }

        if (fmax < ftol) {
            solved = true;
            break;
        }
    
        t_begin = clock();
        S4 a = jac(x);
        t_jac += double(clock() - t_begin) / CLOCKS_PER_SEC;

        NumericVector dims = a.slot("Dim");
        IntegerVector Ap = a.slot("p");
        IntegerVector Ai = a.slot("i");
        NumericVector Ax = a.slot("x");

        // TODO: check that dims[0] == dims[1] = n
    
        // LU factorization
        t_begin = clock();
        void *Symbolic, *Numeric;
        (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), 
                                    REAL(Ax), &Symbolic, null, null) ;
        (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                    Symbolic, &Numeric, null, info) ;
        umfpack_di_free_symbolic (&Symbolic);
        cond = info[UMFPACK_RCOND];
        t_lu += double(clock() - t_begin) / CLOCKS_PER_SEC;

        t_begin = clock();
        (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                 dx, REAL(f), Numeric, null, null);
        t_solve += double(clock() - t_begin) / CLOCKS_PER_SEC;
    
        for (int i = 0; i < n; i++) {
           x[i] -= dx[i];
        }

        umfpack_di_free_numeric(&Numeric);
    }

    if (!control["silent"]) {
        if (solved) {
            Rprintf("Convergence after %d iterations\n", iter);
        } else {
            Rprintf("No convergence after %d iterations\n", iter);
        }
    }

    delete[] dx;

    return List::create(Named("solved")  = solved,
                        Named("iter")    = iter,
                        Named("x")       = x,
                        Named("fval")    = f,
                        Named("t_f")     = t_f,
                        Named("t_jac")   = t_jac,
                        Named("t_lu")    = t_lu,
                        Named("t_solve") = t_solve);
}
