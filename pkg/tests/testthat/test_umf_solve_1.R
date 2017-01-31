library(umfpackr)
library(testthat)

context("Test umf_solve with simple matrix")

# create a sparse matrix
ap <- c(0, 2, 5, 9, 10, 12)
ai <- c(1, 2, 1, 3, 5, 2, 3, 4,  5, 3, 2, 5)
ax <- c(2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1)
b  <- c(8, 45, -3, 3, 19)
a  <- sparseMatrix(i = ai, p = ap,  x = ax,  dims = c(5, 5))

test_that("result of umf_solve is correct", {
    expect_equal(umf_solve(a, b), 1:5)
})
