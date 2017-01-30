// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// umf_solve
NumericVector umf_solve(S4 a, NumericVector b);
RcppExport SEXP umfpackr_umf_solve(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(umf_solve(a, b));
    return rcpp_result_gen;
END_RCPP
}
// umf_solve_nl_
List umf_solve_nl_(NumericVector start, Function fn, Function jac, List control);
RcppExport SEXP umfpackr_umf_solve_nl_(SEXP startSEXP, SEXP fnSEXP, SEXP jacSEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< Function >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< Function >::type jac(jacSEXP);
    Rcpp::traits::input_parameter< List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(umf_solve_nl_(start, fn, jac, control));
    return rcpp_result_gen;
END_RCPP
}