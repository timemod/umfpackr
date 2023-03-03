// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// scale_mat_col
NumericVector scale_mat_col(S4 mat, NumericVector scale);
RcppExport SEXP _umfpackr_scale_mat_col(SEXP matSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(scale_mat_col(mat, scale));
    return rcpp_result_gen;
END_RCPP
}
// umf_solve_
List umf_solve_(S4 a, NumericVector b, List umf_control, bool rowscal);
RcppExport SEXP _umfpackr_umf_solve_(SEXP aSEXP, SEXP bSEXP, SEXP umf_controlSEXP, SEXP rowscalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< List >::type umf_control(umf_controlSEXP);
    Rcpp::traits::input_parameter< bool >::type rowscal(rowscalSEXP);
    rcpp_result_gen = Rcpp::wrap(umf_solve_(a, b, umf_control, rowscal));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_umfpackr_scale_mat_col", (DL_FUNC) &_umfpackr_scale_mat_col, 2},
    {"_umfpackr_umf_solve_", (DL_FUNC) &_umfpackr_umf_solve_, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_umfpackr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
