pure_newton_step <- function(x, dx, fn, ...) {

  x_new <- x + dx
  Fx_new <- fn(x_new, ...)

  return(list(x_new = x_new, Fx_new = Fx_new))
}
