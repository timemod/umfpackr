pure_newton_step <- function(x, dx, iter, cond, fn, control) {

  x_new <- x + dx
  Fx_new <- fn(x_new)

  if (control$trace) {
    report_pure_newton(iter, cond, Fx_new)
  }

  return(list(x_new = x_new, Fx_new = Fx_new))
}


# print an iteration report for pure newton steps
report_pure_newton <- function(iter, cond, Fx) {

  if (iter == 0) {
    cat(sprintf("%5s%15s%20s%20s\n", "Iter", "Jac", "Largest |f|",
                "Index largest |f|"))
  }

  Fx_abs <- abs(Fx)
  Fx_max <- max(Fx_abs)
  if (is.na(Fx_max)) {
    i_max <- Position(is.na, Fx_abs)
  } else {
    i_max  <- which.max(Fx_abs)
  }


  if (iter == 0) {
    cat(sprintf("%5d%15s%20.3e%20d\n", iter, "", Fx_max, i_max))
  } else {
    cat(sprintf("%5d%15.2e%20.3e%20d\n", iter, cond, Fx_max, i_max))
  }
}

