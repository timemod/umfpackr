#include <Rcpp.h>
#include <umfpack.h>

void set_umf_control(double control[UMFPACK_CONTROL], Rcpp::List umf_control);
