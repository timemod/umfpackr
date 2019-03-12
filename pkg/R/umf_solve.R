#' Solves the system of linear equations \eqn{A x = b} using UMFPACK
#'
#' @param a an object of class  \code{dgCMatrix}
#' (see \code{\link[Matrix]{dgCMatrix-class}})
#' @param b the vector \eqn{b}
#' @return the solution \eqn{x}
#' @export
#' @useDynLib umfpackr
#' @importFrom Rcpp sourceCpp
umf_solve <- function(a, b) {
  return(umf_solve_(a, b, 0)$x)
}
