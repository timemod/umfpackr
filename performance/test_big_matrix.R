library(umfpackr)
library(Matrix)

load("example_input.RData")
#print(jac)

cat("\nSolution of umf_solve\n")
print(system.time(
    x_umf <- umf_solve(jac, fval)
))
quit()
cat("\nSolution of Matrix::solve\n")
print(system.time(
    x_matrix <- solve(jac, fval)
))
