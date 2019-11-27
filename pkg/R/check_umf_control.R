check_umf_control <- function(umf_control) {

  if (is.null(umf_control)) {
    return(umf_control)
  } else if (length(umf_control) == 0 && is.list(umf_control)) {
    return(umf_control)
  }

  if (!is.list(umf_control) || is.null(names(umf_control))) {
    stop("Argument umf_control should be a named list.")
  }

  lengths <- sapply(umf_control, FUN = length)
  if (any(lengths == 0 | lengths > 1)) {
    stop("Argument umf_control should be a list of scalar variables.")
  }

  return(umf_control)
}
