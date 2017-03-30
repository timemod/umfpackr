# function to be minimized in the global strategy
get_fnorm <- function(f) {
  return(0.5*sum(f^2))
}
