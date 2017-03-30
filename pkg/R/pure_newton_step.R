pure_newton_step <- function(x, dx, fn, ...) {

  x_new <- x + dx
  f_new <- fn(x_new, ...)

  return(list(x_new = x_new, f_new = f_new))
}
