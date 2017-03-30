library(umfpackr)
library(Matrix)

load("example_input.RData")
#print(jac)
  
print(class(jac))
cat("\ntranspose of jac\n")
print(system.time(
    jac2 <- t(jac)
))

quit()

cat("\nmultiplying\n")
nc <- ncol(jac2)
print(nc)
print(system.time(
    jac3 <- jac %*% numeric(nc)
))


cat("\nSolution of umf_solve\n")
print(system.time(
    x_umf <- umf_solve(jac, fval)
))

cat("\nSolution of Matrix::solve\n")
print(system.time(
    x_matrix <- solve(jac, fval)
))
