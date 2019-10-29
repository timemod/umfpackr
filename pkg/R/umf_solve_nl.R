#' Solves a system of non-linear equations using UMFPACK
#'
#' This function solves a system of non-linear equations  \eqn{F(x) = 0}
#' using Newton's method. UMFPACK is employed to solve the linear equations
#' in each Newton iteration. Optionally a cubic line search is used
#' when a Newton step does not yield a sufficient reduction of the function values.
#' @param start initial guess of the solution \eqn{x}
#' @param fn the function \eqn{F}
#' @param jac a function returning the Jacobian of the function
#' as a \code{dgCMatrix} object
#' @param ... arguments passed to \code{fn} and \code{jac}
#' @param control a list with control parameters. See Details.
#' @param global The global strategy. Possible values are \code{"no"}
#' (no global strategy, the default) and \code{"cline"} (cubic line search)
#' (cubic line search)
#' @param scaling Scaling method. Possible values are
#' \code{"row"}. \code{"col"} and \code{"none"}. The default is \code{"row"}.
#' See Details.
#' @return a list with the following components:
#' \item{\code{solved}}{A logical equal to \code{TRUE} if convergence
#' of the function values has been achieved.}
#' \item{\code{iter}}{the number of iterations}
#' \item{\code{x}}{the final values of \eqn{x}}
#' \item{\code{fval}}{the function value }
#' \item{\code{message}}{A string equal to \code{"ok"} if a solution
#' has been found. Otherwise it describes the reason why the iteration
#' was stopped without success}
#' @details
#' \subsection{Control options}{
#' Argument \code{control} is a named list containing one or more of
#' the following components:
#' \describe{
#' \item{\code{ftol}}{The function value tolerance. Convergence is reached
#' if the largest function value is smaller than \code{ftol}. The default
#' value is \code{1e-8}.}
#' \item{\code{xtol}}{The relative step size tolerance. When the relative
#' step size is smaller than \code{xtol}, then the iteration is stopped.
#' The default value is \code{1e-8}.}
#' \item{\code{maxiter}}{The maximum number of iterations. The default is 20.}
#' \item{\code{trace}}{A logical. If \code{TRUE}  then the progress of the
#' iteraton is printed. The default is \code{FALSE}.}
#' \item{\code{silent}}{A logical. If \code{TRUE}  then all output is suppressed.
#' The default is \code{FALSE}.}
#' \item{\code{allow_singular}}{A logical value (default \code{FALSE})
#' indicating if a small correction to the Jacobian is applied when it is
#' singular or too ill-conditioned.
#' The method used is similar to a Levenberg-Marquardt correction
#' and is explained in Dennis and Schnabel (1996) on page 151.
#' }
#' \item{\code{allow_singular}}{A logical value (default \code{FALSE})
#' indicating if a small correction to the Jacobian is applied when it is
#' singular. The method used is similar to a Levenberg-Marquardt correction
#' and is explained in Dennis and Schnabel (1996) on page 151.
#'}}}
#'\subsection{Scaling of the Jacobian}{
#' For each iteration in the Newton method the linear system \eqn{J s = F(x)} is
#' solved, where the Jacobiab matrix \eqn{J} are the derivatives of the equations
#' with respect to the variables, and \eqn{s} the Newton step.  Scaling
#' can improve the condition of the Jacobian.
#' For \emph{row scaling}, the system is transformed to
#' \eqn{D^{-1} J s = D^{-1} F(x)}, where \eqn{D} is a diagonal matrix with
#' row scaling factors. Here  the scaling factors are the L1 norms of the rows
#' of \eqn{J}. For \emph{column scaling}, the system is transformed to
#' \eqn{J D^{-1} D s = F(x)}, where \eqn{D} is a diagonal matrix with column
#' scaling factors, calculated from  the L1 norms of the columns of \eqn{J}.
#'
#' The scaling is only used to solve the non-linear equations and has no effect
#' on the convergence of the Newton algorihtm. Thus the iterations
#' are considered to be converged when the maximum value of the unscaled
#' function values \eqn{F(x)} is smaller than \code{ftol}.
#' }
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
#' @examples
#'library(umfpackr)
#'
#'dslnex <- function(x, c) {
#'    y <- numeric(2)
#'    y[1] <- x[1]^2 + x[2]^2 - c
#'    y[2] <- exp(x[1]-1) + x[2]^3 - c
#'    y
#'}
#'
#'jacdsln <- function(x, c) {
#'    n <- length(x)
#'    Df <- matrix(numeric(n*n),n,n)
#'    Df[1,1] <- 2*x[1]
#'    Df[1,2] <- 2*x[2]
#'    Df[2,1] <- exp(x[1]-1)
#'    Df[2,2] <- 3*x[2]^2
#'
#'    return(as(Df, "dgCMatrix"))
#'}
#'
#'
#'xstart <- c(2,3)
#'print(umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
#'                   control = list(trace = TRUE)))
#' @importFrom Matrix t
#' @importFrom Matrix norm
#' @importFrom Matrix Diagonal
#' @export
umf_solve_nl <- function(start, fn, jac, ..., control,
                         global = c("no", "cline"),
                         scaling = c("row", "col", "none")) {

  #
  # test arguments
  #
  global <- match.arg(global)
  scaling <- match.arg(scaling)
  if (!is.numeric(start)) {
    stop("Argument 'start' must be a numeric vector.")
  }
  if (!is.function(fn) || !is.function(jac)) {
    stop("Argument 'fn' and 'jac' should be functions.")
  }
  if (!missing(control)) {
    if  (!is.list(control) || (length(control) > 0 &&
                               is.null(names(control)))) {
      stop("Argument 'control' should be a named list.")
    }
  }

  message <- "???"

  control_ <- list(ftol = 1e-8, xtol = 1e-8, maxiter = 20,
                   allow_singular = FALSE,
                   trace = FALSE, silent = FALSE)

  if (!missing(control)) {
      control_[names(control)] <- control
  }

  if (control_$silent) control_$trace <- FALSE

  fun <- function(x) {
    return(fn(x, ...))
  }
  jacfun <- function(x) {
    return(jac(x, ...))
  }

  solved <- FALSE

  if (control_$trace) {
    cat("\nIteration report\n")
    cat("----------------\n")
  }

  n <- length(start)

  x <- start
  cond <- NA_real_
  iter <- 0


  Fx <- fun(x)

  if (!is.numeric(Fx) || length(Fx) != n) {
    stop(sprintf(paste("Function 'fn' should return a numeric vector with the",
                       "same length as argument start (%d).\n"), n))
  }

  colscal <- scaling == "col"
  rowscal <- scaling == "row"

  # initialize scale with zeros
  if (colscal) scale <- numeric(n)

  while (TRUE) {

    Fx_abs <- abs(Fx)
    Fx_max <- max(Fx_abs)

    if (!is.finite(Fx_max)) {
      message <- handle_not_finite_fval(Fx, iter)
      break
    }

    if (iter == 0 && control_$trace) {
      if (global == "cline") {
        report_cline(iter, cond, FALSE, 1, Fx)
      } else if (global == "no") {
        report_pure_newton(iter, cond, Fx)
      }
    }

    if (Fx_max < control_$ftol) {
      solved <- TRUE
      break
    }

    if (iter > 0 && get_step_crit(dx, x) < control_$xtol) {
        solved <- FALSE
        message <- sprintf("Relative step size smaller than xtol (%g)\n",
                           control_$xtol)
        break
    }

    if (iter >= control_$maxiter) {
      message <- sprintf(paste("The maximum number of iterations (%d) has been",
                               "reached\n"), control_$maxiter)
      break
    }

    iter <- iter + 1


    # do new newton step
    j <- jacfun(x)

    if (iter == 1) {
      if (!inherits(j, "dgCMatrix") || nrow(j) != n || ncol(j) != n) {
        stop(sprintf(paste("Function 'jac' should return a dgCMatrix with %d",
                          "rows and columns.\n"), n))
      }
    }

    if (scaling == "col") {
      scale <- scale_mat_col(j, scale)
    }

    # call umf_solve_  to solve j  x = -Fx
    sol <- umf_solve_(j, Fx, rowscal)

    # use a rough estimate of the condition number of UMFPACK.
    cond <- sol$cond

    if (sol$status == "singular matrix") {

      if (any(!is.finite(j@x))) {
        message <- sprintf(paste("The Jacobian contains non-finite values at",
                                 "iteration %d.\n"), iter)
        break
      }

      if (control_$allow_singular) {
        # Use a small perturbation of the Jacobian, See Dennis and Schnabel
        # (1996) on page 151
        n <- length(x)
        h <- t(j) %*% j
        mu <- sqrt(n * .Machine$double.eps) * norm(h, type = "1")
        h <- h + mu * Diagonal(n)
        b <- as.numeric(t(j) %*% Fx)

        sol <- umf_solve_(h, b, rowscal)

        if (sol$status == "singular matrix") {
          message <- sprintf(
                paste("The perturbed Jacobian is still singular.",
                      "The inverse condition is %g.\n"), sol$cond)
          break
        }

     } else {

       message <- sprintf(paste("The Jacobian is singular at iteration %d.",
                                 "The inverse condition is %g.\n"), iter, cond)
       break
     }
    }

    dx <- -sol$x
    if (scaling == "col") dx = dx / scale

    if (global == "no") {
      ret <- pure_newton_step(x, dx, iter, cond, fun, control_)
    } else if (global == "cline") {
      g <- t(j) %*% Fx
      ret <- cline(x, Fx, g, dx, iter, cond, fun, control_)
    }

    if (is.null(ret)) {
      message <- sprintf("No better point found at iter %d\n", iter)
      break
    }

    dx <- x - ret$x_new
    x  <- ret$x_new
    Fx <- ret$Fx_new
  }

  if (solved) {
    message <- "ok"
  }

  if (!control_$silent) {
    if (solved) {
      cat(sprintf("Convergence after %d iterations\n", iter))
    } else {
      cat(message)
    }
  }

  return(list(solved = solved, iter = iter, x = x, fval = Fx,
              message = message))
}


handle_not_finite_fval <- function(f, iter) {
  first_na <- Position(function(x) {x}, !is.finite(f))
  if (iter == 0) {
    return(sprintf(paste("Initial value of function contains",
            "non-finite values (starting at index=%d)\n"), first_na))
  } else {
    return(sprintf(paste("Function value contains",
            "non-finite values (starting at index=%d) at iteration %d\n"),
            first_na, iter))
  }
}
