cublic_linesearch <- function(x, Fx, g, dx, fun, control, ...) {

  ALPHA <- 1e-4

  # TODO: compute LAMBDA_MIN as suggested by Dennis and Schnabel
  LAMBDA_MIN <- 1e-8

  deriv <- as.numeric(t(g) %*% dx)

  lambda <- 1
  lambda_prev <- NA_real_

  f_0 <- get_fnorm(Fx)
  first <- TRUE

  while (TRUE) {

    x_new <- x + lambda * dx
    Fx_new <- fun(x_new, ...)
    # TODO: check NA values
    f_lambda <- get_fnorm(Fx_new)

    if (control$trace) {
      if (first) {
        cat(sprintf(paste0("%20s", paste(rep("#", 80), collapse = ""), "\n"),
                    ""))
        cat(sprintf("%20s Cublic line search\n", ""))
        cat(sprintf("%20s %15s%20s%20s\n", "", "Lambda",
                    "Largest |f|", "Index largest |f|"))
      }
      Fx_abs <- abs(Fx_new)
      Fx_max <- max(Fx_abs)
      i_max <- which.max(Fx_abs)
      cat(sprintf("%20s %15.2e%20.3e%20d\n", "", lambda, Fx_max, i_max))
     }

    if (f_lambda <= f_0 + ALPHA * lambda * deriv) {
      # satisfactory x found
      break;
    } else if (lambda < LAMBDA_MIN)
      # no satisfactory x_new can be found sufficiently distinct from x
      return(NULL)
    else {
      if (first) {
        # first is quadratic
        lambda_new <- - deriv / (2 * (f_lambda - f_0 - deriv))
        first <- FALSE
      } else {
        fac      <- f_lambda       - f_0 - lambda * deriv
        fac_prev <- f_lambda_prev  - f_0 - lambda_prev * deriv
        a <- fac / lambda^2 - fac_prev / lambda_prev^2
        b <- -lambda_prev * fac / lambda^2 + lambda * fac_prev / lambda_prev^2
        cat("a en b bij 1\n")
        print(a)
        print(b)
        a <- a / (lambda - lambda_prev)
        b <- b / (lambda - lambda_prev)
        # TODO: check situation a approx. 0
        disc <- b^2 - 3 * a * deriv
        t1 <- -(b + abs(b) * sqrt(disc)) / (3 * a)
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

    lambda <- max(lambda_new, lambda / 10)
  }

  if (control$trace) {
    cat(sprintf(paste0("%20s", paste(rep("#", 80), collapse = ""), "\n"),
        ""))
  }

  return(list(x_new = x_new, Fx_new = Fx_new))
}
