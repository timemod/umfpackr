# function to be minimized in the global strategy
get_fnorm <- function(f) {
  return(0.5*sum(f^2))
}

# calculate max(abs(dx[*]) / max(x[*], 1))
get_step_crit <- function(dx, x) {
  return(max(abs(dx) / pmax(abs(x), 1)))
}
