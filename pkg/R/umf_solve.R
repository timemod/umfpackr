#' Solves the system of linear equations using UMFPACK
#'
#' This function solved the linear equations of the form \eqn{A x = b} using
#' UMFPACK.
#' @param a an object of class  \code{dgCMatrix}
#' (see \code{\link[Matrix]{dgCMatrix-class}})
#' @param b the vector \eqn{b}
#' @param umf_control  A named list with control parameters passed to UMFPACK.
#' See Details.
#' @return the solution \eqn{x}
#' @references
#' Dennis, J.E. Jr and Schnabel, R.B. (1997), \emph{Numerical Methods for Unconstrained Optimization
#'  and Nonlinear Equations}, Siam.
#'
#' Davis, T.A. (2004). A column pre-ordering strategy for the unsymmetric-pattern
#' multifrontal method. \emph{ACM Trans. Math. Softw.}, \bold{30(2)}, 165–195.
#'
#' Davis, T.A (2004). Algorithm 832: UMFPACK, an unsymmetric-pattern multifrontal
#' method. \emph{ACM Trans. Math. Softw.}, \bold{30(2)}, 196–199.
#'
#' Davis, T.A and Duff, I.S. (1997). An unsymmetric-pattern multifrontal method for
#' sparse LU factorization. \emph{SIAM J. Matrix Anal. Applic.}, \bold{18(1)}, 140–158.
#'
#' Davis, T.A  and Duff, I.S (1999). A combined unifrontal/multifrontal method for
#' unsymmetric sparse matrices. \emph{ACM Trans. Math. Softw.}, \bold{25(1)}, 1–19..
#' @export
#' @useDynLib umfpackr
#' @details
#' \subsection{UMFPACK control options}{
#' Argument `umf_control` can be used to specify UMFPACK control parameters.
#' Consult the documentation of UMFPACK for a description of the various control
#' paramters. `umf_control` should be a named list; the names are the
#' names of the UMFPACK control parameters excluding the prefix `UMFPACK`.
#' With a few exceptions described below, the values are numerical values.
#' For example,
#' ```
#' list(SYM_PIVOT_TOLERANCE = 0.01, AGGRESIVE = 1)
#' ```
#' The values for the control options `STRATEGY`, `ORDERING` and `SCALE`
#' should be character. In UMFPACK the values of these parameters
#' should be specified with named numerical constants. For example,
#' for `UMFPACK_STRATEGY` allowed values are the constants
#'  `UMFPACK_STRATEGY_AUTO`, `UMFPACK_STRATEGY_UNSYMMETRIC` and
#'  `UMFPACK_STRATEGY_SYMMETRIC`. In package `umfpackr` the values should be
#'  the name of the corresponding constant , again excluding the prefix `UMFACK_`.
#' Exanple:
#' ```
#' list(STRATEGY = "STRATEGY_UNSYMMETRIC",
#'      ORDERING = "ORDERING_METIS",
#'      SCALE    = "SCALE_NONE")
#' ```
#' }
#' @importFrom Rcpp sourceCpp
umf_solve <- function(a, b, umf_control = list()) {

  if (!inherits(a, "dgCMatrix")) {
    stop("a is not an object of class dgCMatrix")
  }

  if (!is.numeric(b)) {
    stop("b is not a numerical vector or matrix")
  }

  if (is.matrix(b) && ncol(b) > 1) {
    stop("b should be a numeric vector of a matrix with 1 column")
  }

  if (length(b) != nrow(a)) {
    stop(sprintf(paste("The length of vector b (%d) is not equal to the",
                       "number of rows of a (%d)."), length(b), nrow(a)))
  }

  umf_control <- check_umf_control(umf_control)

  sol <- umf_solve_(a, b, umf_control)
  if (sol$status == "singular matrix") {
    stop(sol$status)
  }
  return(sol$x)
}
