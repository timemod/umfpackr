#' Solves the system of linear equations \eqn{A x = b} using UMFPACK
#'
#' @param a a \code{\link[Matrix]{dgCMatrix}} 
#' @param b the vector \eqn{b}
#' @return the solution \eqn{x}
#' @export
#' @useDynLib umfpackr
#' @importFrom Rcpp sourceCpp
umf_solve <- function(a, b) {
  return(umf_solve_(a, b)$x)
}
