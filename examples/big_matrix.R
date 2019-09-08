library(umfpackr)
library(Matrix)

load("input/example_input.RData")
#print(jac)

cat("\nSolution of umf_solve\n")
print(system.time(
  x_umf <- umf_solve(jac, fval)
))

cat("\nSolution of Matrix::solve\n")
print(system.time(
  x_matrix <- solve(jac, fval)
))
