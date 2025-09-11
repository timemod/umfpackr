# This function reads the  inverse rcondition of the Jacobian matrix
# from the output of umpack_solve_nl
read_rconds <- function(report) {
  itr_lines <- report[6:(length(report) - 1)]
  ma <- stringr::str_match(itr_lines, "^\\s*\\d+\\s+(.+?)\\s")
  rconds <- as.numeric(ma[, 2])
  return(rconds)
}
