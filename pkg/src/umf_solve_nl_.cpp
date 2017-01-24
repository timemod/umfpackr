#include <Rcpp.h>
#include "umfpack.h"
using namespace Rcpp;

// [[Rcpp::export]]
List umf_solve_nl_(NumericVector start, Function fn, Function jac,
                   List control) {

    bool trace    = control["trace"];
    double ftol   = control["ftol"];
    int maxiter   = control["maxiter"];
    double cnmtrx = control["cnmtrx"];

    //initialize initial convergence speed threshold

    bool solved = false;

    if (trace) {
        Rprintf("\nIteration report\n");
        Rprintf("----------------\n");
        Rprintf("%5s%20s%20s%20s%20s\n", "Iter", "inv. cond. jac.",
                    "Largest |f|", "Index largest |f|", "Fcrit");
    }

    double cond;
    int iter;

    // fmax  largest scaled change in feedback variables
    // fquot  ratio of current and previous fmax
    // fcrit  geometric mean current and previous Fquot
    double fstart = cnmtrx > 0.5 ?  0.5 : 0.5 * cnmtrx;
    double fquot  = fstart;
    double fquotp = fstart;
    double fcrit  = fstart;
    double fmax = 0, fmaxp = 0;
    bool calc_matrix = true;
    bool matrix_calculated = false;

    void *Symbolic = NULL, *Numeric = NULL;
    IntegerVector dims, Ap, Ai;
    NumericVector Ax;
    double *null = (double *) NULL;

    int n = start.size();
    NumericVector f;
    NumericVector x(n), dx(n);
    x = start;


    for (iter = 0; iter < maxiter; iter++) {

        fmaxp = fmax;
        fquotp = fquot;

        f = fn(x);
        NumericVector f_abs = abs(f);
        fmax = max(f_abs);
        int i = which_max(f_abs);
        // TODO: handle NA values

        if (matrix_calculated) {
            Rprintf("%5d%20.2e%20.3e%20d%20.2f\n", iter, cond, fmax, i, fcrit);
        } else {
            Rprintf("%5d%20s%20.3e%20d%20.2f\n", iter, "", fmax, i, fcrit);
        }

        if (fmax < ftol) {
            solved = true;
            break;
        }

        if (iter > 0) {
            fquot = fmax / fmaxp;
            fcrit = sqrt(fquot * fquotp);
            calc_matrix = fcrit > cnmtrx;
        }


        if (calc_matrix) {

            // TODO: separate function
    
            if (Numeric != NULL) {
                umfpack_di_free_numeric (&Numeric) ;
            }
            S4 a = jac(start);
            dims = a.slot("Dim");
            Ap = a.slot("p");
            Ai = a.slot("i");
            Ax = a.slot("x");
    
            // TODO: check that dims[0] == dims[1] = n?
    
            // LU factorization
            (void) umfpack_di_symbolic (n, n, INTEGER(Ap), INTEGER(Ai), 
                                        REAL(Ax), &Symbolic, null, null) ;
            (void) umfpack_di_numeric (INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                       Symbolic, &Numeric, null, null) ;
            umfpack_di_free_symbolic (&Symbolic) ;
    
            // TODO: compute the inverse condition of the matrix
            cond = NA_REAL;
        
            calc_matrix = false;
            matrix_calculated = true;
        } else {
            matrix_calculated = false;
        }

        // solve
        (void) umfpack_di_solve (UMFPACK_A, INTEGER(Ap), INTEGER(Ai), REAL(Ax),
                                 REAL(dx), REAL(f), Numeric, null, null) ;

        x = x - dx;
    }

    if (Numeric != NULL) {
        umfpack_di_free_numeric(&Numeric);
    }

    return List::create(Named("solved") = solved,
                        Named("iter")   = iter,
                        Named("x")      = x,
                        Named("fval")   = f);
}
