library(umfpackr)
library(testthat)

rm(list = ls())

context("Test umf_solve with simple matrix")

# create a sparse matrix
ap <- c(0, 2, 5, 9, 10, 12)
ai <- c(1, 2, 1, 3, 5, 2, 3, 4,  5, 3, 2, 5)
ax <- c(2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1)
b  <- c(8, 45, -3, 3, 19)
a  <- sparseMatrix(i = ai, p = ap,  x = ax,  dims = c(5, 5))

test_that("result of umf_solve is correct", {
  expect_equal(umf_solve(a, b), 1:5)
  expect_equal(umf_solve(a, as.matrix(b)), 1:5)
})

test_that("singular matrix", {
  a_singular <- a
  a_singular[ , 1] <- 0
  expect_error(umf_solve(a_singular, b), "singular matrix")
})

test_that("several errors", {
  expect_error(umf_solve(b, b), "a is not an object of class dgCMatrix")
  expect_error(umf_solve(a, cbind(b, b)),
          "b should be a numeric vector of a matrix with 1 column")
  expect_error(umf_solve(a, b[1:2]),
        "The length of vector b \\(2\\) is not equal to the number of rows of a \\(5\\).")
})
