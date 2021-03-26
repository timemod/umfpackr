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
#' (no global strategy, the default) and \code{"cline"} (cubic line search).
#' @param scaling Scaling method. Possible values are
#' \code{"row"}. \code{"col"} and \code{"none"}. The default is \code{"row"}.
#' See Details.
#' @param umf_control A named list with control parameters passed to UMFPACK.
#' Currently only a single control parameter can specified: `ordering`, which
#' specifies the ordering method. Allowed values are `"AMD"` (the default),
#' `"CHOLMOD"`, `"METIS"` and `"BEST"`. See the Vignette
#' \href{../doc/UMFPACK_interface.pdf}{\emph{UMFPACK interface for R}} for more
#' details. For example, to use METIS ordering specify
#' `umf_control = list(ORDERING = "METIS")`. The METIS ordering method can
#' handle larger matrices than the standard AMD method.
#' @return a list with the following components:
#' \item{\code{solved}}{A logical equal to \code{TRUE} if convergence
#' of the function values has been achieved.}
#' \item{\code{iter}}{the number of iterations}
#' \item{\code{x}}{the final values of \eqn{x}}
#' \item{\code{fval}}{the final  function values}
#' \item{\code{message}}{A string equal to \code{"ok"} if a solution
#' has been found. Otherwise it describes the reason why the iteration
#' was stopped without success}
#' @details
#' \subsection{Control options}{
#' Argument \code{control} is a named list containing one or more of
#' the following components:
#' \describe{
#' \item{\code{ftol}}{The function value tolerance. Convergence is reached
#' if the largest absolute function value is smaller than \code{ftol}.
#' The default value is \code{1e-8}.}
#' \item{\code{xtol}}{The relative step size tolerance. If convergence
#' has not been reached yet and if the relative
#' step size of all \code{x} values is smaller than \code{xtol} then the
#' iteration process is terminated with an error.
#' The default value is \code{1e-8}. The relative step size of
#' \eqn{x[i]} is calculated as
#' \eqn{|x[i] - x^*[i]| / \rm{max}(|x[i]|, 1)}, where \eqn{x^*[i]} is
#' the value of \eqn{x[i]} at the previous iteration.}
#' \item{\code{maxiter}}{The maximum number of iterations. The default is 20
#' if no global strategy is used (argument `global = "no"`), and 150
#' if cubic line searching is used (argument `global = "cline"`).}
#' \item{\code{trace}}{A logical. If \code{TRUE}  then the progress of the
#' iteration is printed. The default is \code{FALSE}.}
#' \item{\code{silent}}{A logical. If \code{TRUE}  then all output is suppressed.
#' The default is \code{FALSE}.}
#' \item{\code{cnd_tol}}{The tolerance for the inverse condition of the jacobian.
#' If the inverse condition is smaller than `cnd_tol`, the solution process
#' is terminated with an error, except if control parameter `allow_singular` is
#' set to `TRUE`. The default is the machine precision (on most platforms
#' about `2e-16`). If the inverse condition is very small but nonzero
#' it may be difficult to find a solution, or the solution may not be
#' meaningful. However, sometimes a good solution can be found even if the
#' condition is quite small. The test  for the ill-conditioning of the
#' jacobian can be turned off by setting `cnd_tol` to 0 or a negative number.
#' However, if the matrix is singular (the inverse condition is exactly zero),
#' it is never possible to continue with the solution process
#'  (except if control parameter `allow_singular` is set to `TRUE`).
#'  The default value of `cnd_tol` is quite small, in some cases it can be
#'  appropriate to use a somewhat larger value (for example `1e-12`)}
#' \item{\code{cnd_method}}{A character vector specifying the method used to
#' estimate the inverse condition number of the jacobian. Possible options are
#' `"umfpack"`(the default), `"condest"` and `"kappa"`.
#' For `"umfpack"` a rough estimate of the condition as computed by UMFPACK is
#' used, using the expression
#' \eqn{\rm{min}(\rm{abs}(\rm{diag}(U)))/\rm{max}(\rm{abs}(\rm{diag}(U)))},
#' where \eqn{U} is the \eqn{U} matrix of the LU factorisation of the jacobian.
#' Method `"condest"` employs function \code{\link[Matrix]{condest}} of the `Matrix`
#' package and `kappa` the function \code{\link[base]{kappa}} of the `base` package.
#' Method `condtest` is more accurate than the rough estimate of UMFPACK,
#' but takes more time. `kappa` is exact, but is slow for large
#' matrices because this function does not use sparse matrices.
#' Method `"condest"` usually gives a reasonable approximation of the
#' inverse condition number.
#' It is recommended to normally use `"umfpack"`, but occasionally use `"condest"`
#' for a more accurate check of the condition number.}
#' \item{\code{allow_singular}}{A logical value (default \code{FALSE})
#' indicating if a small correction to the Jacobian is applied when it is
#' singular or too ill-conditioned.
#' The method used is similar to a Levenberg-Marquardt correction
#' and is explained in Dennis and Schnabel (1996) on page 151.
#' If `TRUE`, then the correction is applied if the inverse condition is
#' exactly zero or if the inverse condition is smaller than control
#' parameter `cnd_tol`.
#' }
#' }}
#'\subsection{Scaling of the Jacobian}{
#' For each iteration in the Newton method the linear system \eqn{J s = F(x)} is
#' solved, where the Jacobian matrix \eqn{J} are the derivatives of the equations
#' with respect to the variables, and \eqn{s} the Newton step.  Scaling
#' can improve the condition of the Jacobian.
#' For \emph{row scaling}, the system is transformed to
#' \eqn{D^{-1} J s = D^{-1} F(x)}, where \eqn{D} is a diagonal matrix with
#' row scaling factors. Here  the scaling factors are the L1 norms of the rows
#' of \eqn{J}. For \emph{column scaling}, the system is transformed to
#' \eqn{J D^{-1} D s = F(x)}, where \eqn{D} is a diagonal matrix with column
#' scaling factors, calculated from  the L1 norms of the columns of \eqn{J}.
#'
#' The scaling is only used to solve the linear equations and has no effect
#' on the convergence of the Newton algorithm. Thus the iterations
#' are considered to be converged when the maximum value of the unscaled
#' function values \eqn{F(x)} is smaller than \code{ftol}.
#' It is therefore recommended to define a the function \eqn{F(x)} for which all
#' function values have comparable orders of magnitude and all
#' function arguments \eqn{x} have similar orders of magnitude.
#' }
#' @references
#' Dennis, J.E. Jr and Schnabel, R.B. (1997), \emph{Numerical Methods for
#' Unconstrained Optimisation and Nonlinear Equations}, Siam.
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
#'
#' # now use METIs columns ordering (run this on Linux only)
#' if (.Platform$OS.type != "windows") {
#'    print(umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
#'                   control = list(trace = TRUE),
#'                   umf_control = list(ordering = "METIS")))
#' }
#' @importFrom Matrix t
#' @importFrom Matrix norm
#' @importFrom Matrix Diagonal
#' @importFrom Matrix condest
#' @seealso \code{\link{umf_solve}}.
#' @export
umf_solve_nl <- function(start, fn, jac, ..., control,
                         global = c("no", "cline"),
                         scaling = c("row", "col", "none"),
                         umf_control = list()) {
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

  # create default control options:
  maxiter <- if (global == "no") 20 else 150
  control_ <- list(ftol = 1e-8, xtol = 1e-8, maxiter = maxiter,
                   allow_singular = FALSE,
                   cnd_tol = .Machine$double.eps,
                   cnd_method = "umfpack",
                   trace = FALSE, silent = FALSE)

  if (!missing(control)) {
    if  (!is.list(control) || (length(control) > 0 &&
                               is.null(names(control)))) {
      stop("Argument 'control' should be a named list.")
    }
    if (length(unknown_control_options <- setdiff(names(control),
                                                  names(control_))) > 0) {
      stop("Unknown control options ", paste(unknown_control_options,
                                             sep = ", "), ".")
    }
    if (!is.null(control$cnd_method)) {
      if (!is.character(control$cnd_method) ||
          length(control$cnd_method) != 1) {
        stop("Control option 'cnd_method' should be a character of length 1")
      }
      if (!control$cnd_method %in% c("umfpack", "condest", "kappa")) {
        stop("Allowed values for control option 'cnd_method' are",
             " 'umfpack', 'condest' and 'kappa'.")
      }
    }
    if (!is.null(control$cnd_tol)) {
      if (!is.numeric(control$cnd_tol) || length(control$cnd_tol) != 1) {
        stop("Control option 'cnd_tol' should be a numeric of length 1")
      }
    }
    control_[names(control)] <- control
    if (control_$silent) control_$trace <- FALSE
  }

  umf_control <- check_umf_control(umf_control)

  colscal <- scaling == "col"
  rowscal <- scaling == "row"

  fun <- function(x) {
    return(fn(x, ...))
  }
  jacfun <- function(x) {
    return(jac(x, ...))
  }

  # Function for a more accurate calculation of the inverse condition number
  # of the jacobian than the rough estimate provided by UMFPACK.
  get_cond <- function(jac) {
    if (rowscal) {
      # Row scaling is applied internally in UMFPACK. If we want to estimate the
      # condition number with the 'condest' or 'kappa' method, we should
      # first apply row scaling here again.
      scale <- 1 / Matrix::rowSums(abs(jac))
      scale <- ifelse(is.finite(scale), scale, 1)
      jac <- jac * scale
    }
    # Note: if colscal == TRUE, then jac is already scaled, so scaling
    # is not needed.
    if (control_$cnd_method == "condest") {
      cond <- 1 / condest(jac)$est
    } else if (control_$cnd_method == "kappa") {
      cond <- 1 / kappa(jac, exact = TRUE)
    }
    return(cond)
  }

  message <- "???"
  solved <- FALSE

  n <- length(start)

  x <- start
  cond <- NA_real_
  iter <- 0

  Fx <- fun(x)

  if (!is.numeric(Fx) || length(Fx) != n) {
    stop(sprintf(paste("Function 'fn' should return a numeric vector with the",
                       "same length as argument start (%d).\n"), n))
  }

  # Initialize scale with zeros:
  if (colscal) scale <- numeric(n)

  while (TRUE) {

    Fx_abs <- abs(Fx)
    Fx_max <- max(Fx_abs)

    if (!is.finite(Fx_max)) {
      message <- handle_not_finite_fval(Fx, iter)
      break
    }

    if (iter == 0 && control_$trace) {
      cat("\nIteration report\n")
      cat("----------------\n")
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

    # Compute the jacobian:
    j <- jacfun(x)

    if (iter == 1) {
      if (!inherits(j, "dgCMatrix") || nrow(j) != n || ncol(j) != n) {
        stop(sprintf(paste("Function 'jac' should return a dgCMatrix with %d",
                          "rows and columns.\n"), n))
      }
    }

    # Calculate the column scale factors and scale the matrix using C++ function
    # scale_mat_col. Note: objecy j is modified in function scale_mat_col.
    if (colscal) scale <- scale_mat_col(j, scale)
    # Call C++ function umf_solve_ to solve j  x = Fx
    sol <- umf_solve_(j, Fx, umf_control, rowscal)

    # if cnd_method == "umfpack", use a rough estimate of the inverse condition
    # available in the UMFPACK result, otherwise use more advanced methods
    # defined in function get_cond.
    cond <- if (control_$cnd_method == "umfpack") sol$cond else get_cond(j)

    if (sol$status == "singular matrix" || cond < control_$cnd_tol ) {

      if (sol$status == "singular matrix" && any(!is.finite(j@x))) {
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

        sol <- umf_solve_(h, b, umf_control, rowscal)
        cond <- if (control_$cnd_method == "umfpack") sol$cond else get_cond(h)
        if (sol$status == "singular matrix" || cond < control_$cnd_tol) {
          message <- sprintf(paste("The perturbed Jacobian is still singular.",
                                   "The inverse condition is %g.\n"), cond)
          break
        }
      } else {
        if (sol$status == "singular matrix") {
          message <- sprintf(paste("The Jacobian is singular at iteration %d.",
                                   "The inverse condition is %g.\n"),
                             iter, cond)
        } else {
          message <- sprintf(paste("The inverse condition of the jacobian is",
                                   "smaller than cnd_tol (%.3g) at iteration",
                                   "%d.\nThe inverse condition is %.3g.\n"),
                             control_$cnd_tol, iter, cond)
        }
        break
      }
    }

    dx <- -sol$x
    if (colscal) dx = dx / scale

    if (global == "no") {
      ret <- pure_newton_step(x, dx, iter, cond, fun, control_)
    } else if (global == "cline") {
      g <- t(j) %*% Fx
      ret <- cline(x, Fx, g, dx, iter, cond, fun, control_)
      if (is.null(ret)) {
        message <- sprintf("No better point found at iter %d\n", iter)
        break
      }
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
      cat("\n")
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
