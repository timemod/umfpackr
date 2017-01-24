library(umfpackr)
library(Matrix)

ap <- c(0, 2, 5, 9, 10, 12)
ai <- c(1, 2, 1, 3, 5, 2, 3, 4,  5, 3, 2, 5)
ax <- c(2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1)
b <- c(8, 45, -3, 3, 19)


a <- sparseMatrix(i = ai, p = ap,  x = ax,  dims = c(5,5))
print(a)

cat("\nSolution of umf_solve\n")
print(umf_solve(a, b))
cat("\nSolution of Matrix::solve\n")
print(solve(a, b))
