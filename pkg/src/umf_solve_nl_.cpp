#include <Rcpp.h>
#include "umfpack.h"
using namespace Rcpp;

//#define TIMER
#ifdef TIMER
#include <time.h>
clock_t begin, end, t_f_begin, t_f_end, t_jac_begin, t_jac_end,
        t_lu_begin, t_lu_end, t_solve_begin, t_solve_end;
double t_f, t_jac, t_lu, t_solve, delta_t;
#endif

// [[Rcpp::export]]
List umf_solve_nl_(NumericVector start, Function fn, Function jac,
                   List control) {

#ifdef TIMER
    begin = clock();
#endif

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
    double cond;

    int n = start.size();
    NumericVector f;
    NumericVector x(clone(start)), dx(n);

#ifdef TIMER
    t_f = 0.0;
    t_jac = 0.0;
    t_lu = 0.0;
    t_solve = 0.0;
#endif

    for (iter = 0; iter < maxiter; iter++) {

#ifdef TIMER
        t_f_begin = clock();
#endif
        f = fn(x);
#ifdef TIMER
        t_f_end = clock();
        t_f += double(t_f_end - t_f_begin) / CLOCKS_PER_SEC;
#endif
        NumericVector f_abs = abs(f);
        double fmax = max(f_abs);

        // TODO: handle NA values
        if (R_IsNA(fmax)) {
            int i_na;
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
    
#ifdef TIMER
        t_jac_begin = clock();
#endif
        S4 a = jac(x);
        NumericVector dims = a.slot("Dim");
        IntegerVector Ap = a.slot("p");
        IntegerVector Ai = a.slot("i");
        NumericVector Ax = a.slot("x");
#ifdef TIMER
        t_jac_end = clock();
        delta_t = double(t_jac_end - t_jac_begin) / CLOCKS_PER_SEC;
        t_jac +=  delta_t;
#endif

        // TODO: check that dims[0] == dims[1] = n
    
        // LU factorization
#ifdef TIMER
        t_lu_begin = clock();
#endif
        void *Symbolic, *Numeric;
        (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), 
                                    REAL(Ax), &Symbolic, null, null) ;
        (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                    Symbolic, &Numeric, null, info) ;
        umfpack_di_free_symbolic (&Symbolic);
        cond = info[UMFPACK_RCOND];
#ifdef TIMER
        t_lu_end = clock();
        t_lu +=  double(t_lu_end - t_lu_begin) / CLOCKS_PER_SEC;
#endif

#ifdef TIMER
        t_solve_begin = clock();
#endif
        (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                 REAL(dx), REAL(f), Numeric, null, null);
#ifdef TIMER
        t_solve_end = clock();
        t_solve +=  double(t_solve_end - t_solve_begin) / CLOCKS_PER_SEC;
#endif

        x = x - dx;

        umfpack_di_free_numeric(&Numeric);
    }

    if (!control["silent"]) {
        if (solved) {
            Rprintf("Convergence after %d iterations\n", iter);
        } else {
            Rprintf("No convergence after %d iterations\n", iter);
        }
    }

#ifdef TIMER
    end = clock();
    Rprintf("Total time                   %g\n", 
                        double(end - begin) / CLOCKS_PER_SEC);
    Rprintf("timing functon evaluation    %g\n", t_f);
    Rprintf("timing Jacobian calculation: %g\n", t_jac);
    Rprintf("timing LU:                   %g\n", t_lu);
    Rprintf("timing solve:                %g\n", t_solve);

#endif

    return List::create(Named("solved") = solved,
                        Named("iter")   = iter,
                        Named("x")      = x,
                        Named("fval")   = f);
}
