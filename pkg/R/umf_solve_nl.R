#' Solves a system of non-linear equations \eqn{F(x) = 0} using UMFPACK
#'
#' @param start initial guess of the solution \eqn{x}
#' @param fn the function \eqn{F}
#' @param jac a function returning the Jacobian of the function
#' as a \code{dgCMatrix} object
#' @param ... arguments passed to \code{fn} and \code{jac}
#' @param control a list with control parameters
#' @return a list with information about the solution
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
#' @export
umf_solve_nl <- function(start, fn, jac, ..., control = list(),
                         global = c("no", "cline")) {

  global <- match.arg(global)

  control_ <- list(ftol = 1e-8, maxiter = 20, trace = FALSE,
                   cndtol = 1e-12, silent = FALSE)

  control_[names(control)] <- control

  fun <- function(x) {
    return(fn(x, ...))
  }
  jacfun <- function(x) {
    return(jac(x, ...))
  }

  solved <- FALSE

  if (control_$trace) {
    cat("\nIteration report\n");
    cat("----------------\n");
    cat(sprintf("%5s%15s%20s%20s\n", "Iter", "Jac",
            "Largest |f|", "Index largest |f|"))
  }

  x <- start
  cond <- NA_real_
  iter <- 0


  Fx <- fun(x)

  while (TRUE) {

    Fx_abs <- abs(Fx)
    Fx_max <- max(Fx_abs)

    if (is.na(Fx_max)) {
      handle_na_in_fval(Fx, iter)
      break
    }

    if (control_$trace) {
      i <- which.max(Fx_abs);
      if (iter > 1) {
        cat(sprintf("%5d%15.2e%20.3e%20d\n", iter, cond, Fx_max, i))
      } else {
        cat(sprintf("%5d%15s%20.3e%20d\n", iter, "", Fx_max, i))
      }
    }

    if (Fx_max < control_$ftol) {
      solved <- TRUE
      break
    }

    iter <- iter + 1
    if (iter > control_$maxiter) {
      break
    }

    # do new newton step
    j <- jacfun(x)
    sol <-umf_solve_(j, Fx)
    dx <- - sol$x
    cond <- sol$cond

    if (global == "no") {
      ret <- pure_newton_step(x, dx, fun)
    } else if (global == "cline") {
      g <- t(j) %*% Fx
      ret <- cublic_linesearch(x, Fx, g, dx, fun, control_)
    }

    x <- ret$x_new
    Fx <- ret$Fx_new
  }

  if (!control_$silent) {
    if (solved) {
      cat(sprintf("Convergence after %d iterations\n", iter))
    } else {
      cat(sprintf("No convergence after %d iterations\n", iter))
    }
  }

  return(list(solved = solved, iter = iter, x = x, fval = Fx))
}


handle_na_in_fval <- function(f, iter) {
  first_na <- Position(function(x) {x}, is.na(f))
  if (iter == 0) {
    cat(sprintf(paste("Initial value of function contains",
            "non-finite values (starting at index=%d)\n"), first_na))
  } else {
    cat(sprintf(paste("Function value contains",
            "non-finite values (starting at index=%d)\n"), first_na))
  }
  return(invisible(NULL))
}
