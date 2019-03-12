#' Solves a system of non-linear equations \eqn{F(x) = 0} using UMFPACK
#'
#' @param start initial guess of the solution \eqn{x}
#' @param fn the function \eqn{F}
#' @param jac a function returning the Jacobian of the function
#' as a \code{dgCMatrix} object
#' @param ... arguments passed to \code{fn} and \code{jac}
#' @param control a list with control parameters
#' @param global The global strategy. Possible values are \code{"no"}
#' (no global strategy, the default) and \code{"cline"} (cubic line search)
#' (cubic line search)
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
#' @importFrom Matrix norm
#' @importFrom Matrix Diagonal
#' @importFrom Matrix condest
#' @export
umf_solve_nl <- function(start, fn, jac, ..., control = list(),
                         global = c("no", "cline")) {

  global <- match.arg(global)

  message <- "???"

  control_ <- list(ftol = 1e-8, xtol = 1e-8, maxiter = 20,
                   allow_singular = FALSE, trace = FALSE, cndtol = 1e-12,
                   acc_cnd = FALSE, silent = FALSE)

  control_[names(control)] <- control

  control_$cndtol <- max(control_$cndtol, .Machine$double.eps)

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

  x <- start
  cond <- NA_real_
  iter <- 0


  Fx <- fun(x)

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
        solved <- TRUE
        break
    }

    if (iter >= control_$maxiter) {
      message <- sprintf(paste("The maximum number of iterations (%d) has been",
                               "reached"), control_$maxiter)
      break
    }

    iter <- iter + 1


    # do new newton step
    j <- jacfun(x)

    if (control_$acc_cnd) {
      # Estimate the condition number using the Matrix package
      # This may cost considerable CPU time.
      cond <- 1 / condest(j)$est
    } else {
      # call umf_solve_ first to obtain a rough estimate of the condition number.
      # actually, cond < cndtol, we do not accept the solution anyway so we
      # could skip solving the model.
      sol <- umf_solve_(j, Fx, control_$cndtol)
      # use a rough estimate of the condition number of UMFPACK.
      cond <- sol$cond
    }

    if (cond < control_$cndtol) {
      if (!control_$allow_singular) {
        message <- sprintf(paste("The Jacobian is (nearly) singular at",
                          "iteration %d.",
                          "The inverse condition is %g.\n"), iter, cond)
        break
      }
      # Use a small perturbation of the Jacobian, See Dennis and Schnabel (1996)
      # on page 151
      n <- length(x)
      h <- t(j) %*% j
      mu <- sqrt(n * .Machine$double.eps) * norm(h, type = "1")
      h <- h + mu * Diagonal(n)
      b <- as.numeric(t(j) %*% Fx)
      sol <- umf_solve_(h, b, 0)
      if (sol$cond < .Machine$double.eps) {
        # this situation should theoretically not happen
        message <- sprintf(
              paste("The perturbed Jacobian is still (nearly) singular.",
                    "The inverse condition is %g.\n"), sol$cond)
        break
      }

    } else {
      # for acc_cnd, wew still need to solve
      sol <- umf_solve_(j, Fx, 0)
    }

    dx <- - sol$x

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
            "non-finite values (starting at index=%d)\n"), first_na))
  }
}
