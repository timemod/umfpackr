# cubic linesearch
cline <- function(x, Fx, g, dx, iter, cond, fun, control) {

  ALPHA <- 1e-4
  LAMBDA_MIN <- control$xtol / get_step_crit(dx, x)

  deriv <- as.numeric(t(g) %*% dx)

  lambda <- 1
  lambda_prev <- NA_real_

  f_0 <- get_fnorm(Fx)
  first <- TRUE

  while (TRUE) {

    x_new <- x + lambda * dx
    Fx_new <- fun(x_new)
    f_lambda <- get_fnorm(Fx_new)

    if (control$trace) {
        report_cline(iter, cond, first, lambda, Fx_new)
    }

    f_lambda_na <- is.na(f_lambda)

    if (f_lambda_na) {
      lambda_new <- lambda / 2
      first <- FALSE
    } else if (f_lambda <= f_0 + ALPHA * lambda * deriv) {
      # satisfactory x found
      break
    } else if (lambda <= LAMBDA_MIN) {
      # no satisfactory x_new can be found sufficiently distinct from x
      return(NULL)
    } else {
      if (first || f_lambda_prev_na) {
        # first is quadratic
        lambda_new <- - deriv / (2 * (f_lambda - f_0 - deriv))
        if (!first) lambda_new <- min(lambda_new, lambda / 2)
        first <- FALSE
      } else {
        fac      <- f_lambda       - f_0 - lambda * deriv
        fac_prev <- f_lambda_prev  - f_0 - lambda_prev * deriv
        a <- fac / lambda^2 - fac_prev / lambda_prev^2

        b <- -lambda_prev * fac / lambda^2 + lambda * fac_prev / lambda_prev^2
        a <- a / (lambda - lambda_prev)
        b <- b / (lambda - lambda_prev)
        # TODO: check situation a approx. 0

        disc <- b^2 - 3 * a * deriv
        t1 <- -(b + sign(b) * sqrt(disc)) / (3 * a)
        t2 <- deriv / (3 * a) / t1
        if (a > 0) {
          # upward opening parabola ==> rightmost is solution
          lambda_new <- max(t1, t2)
        } else {
          # downward opening parabola ==> leftmost is solution
          lambda_new <- min(t1, t2)
        }
        lambda_new <- min(lambda_new, lambda / 2)
      }
    }

    lambda_prev   <- lambda
    f_lambda_prev <- f_lambda
    f_lambda_prev_na <- f_lambda_na

    lambda <- max(lambda_new, lambda / 10)
  }

  return(list(x_new = x_new, Fx_new = Fx_new))
}

# print an iteration report for the cubic line search
report_cline <- function(iter, cond, jac, lambda, Fx) {

  if (iter == 0) {
    cat(sprintf("%5s%15s%15s%20s%20s\n", "Iter", "Jac", "Lambda", "Largest |f|",
                "Index largest |f|"))
  }

  Fx_abs <- abs(Fx)
  Fx_max <- max(Fx_abs)
  i_max  <- which.max(Fx_abs)

  if (jac) {
    cat(sprintf("%5d%15.2e%15.2e%20.3e%20d\n", iter, cond, lambda, Fx_max,
                i_max))
  } else {
    cat(sprintf("%5d%15s%15.2e%20.3e%20d\n", iter, "", lambda, Fx_max, i_max))
  }
}
