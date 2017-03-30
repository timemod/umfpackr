cublic_linesearch <- function(x, f, g, dx, fn, ...) {

  ALPHA <- 1e-4

  cat("g\n")
  print(g)

  cat("dx\n")
  print(dx)

  slope <- as.numeric(t(g) %*% dx)

  cat("slope\n")
  print(slope)

  lambda <- 1

  fnorm <- get_fnorm(f)
  first <- TRUE

  while (TRUE) {

    x_new <- x + lambda * dx
    f_new <- fn(x_new, ...)
    fnorm_new <- get_fnorm(f_new)

    if (fnorm_new <= fnorm + ALPHA * lambda *slope) {
      break;
    }

    if (first) {
      # first is quadratic
      t <- ((- lambda^2 * slope / 2) / (fnorm_new - fnorm - lambda - slope))
      first <- FALSE
    } else {

    }

    lambda0 <- lambda
    fnorm0 <- fnorm_new
    lambda <- max(t, lambda / 10)


  }

  return(list(x_new = x_new, f_new = f_new))
}
