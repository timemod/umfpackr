library(testthat)
library(umfpackr)
library(Matrix)

context("test for singular matrix (2)")

rm(list = ls())

f <- function(x) {
  return(sqrt(x) + 0.1)
}

jac <- function(x) {
  j <- matrix(1/sqrt(x), nrow = 1, ncol = 1)
  return(as(j, "dgCMatrix"))
}



test_that("starting value 0", {
  expect_known_output(
    ret <- umf_solve_nl(0, f, jac, control = list(trace = TRUE)),
    "expected_output/singular2_report1.txt")
  expect_equal(ret, list(solved = FALSE, iter = 1, x = 0, fval = 0.1,
                         message = paste("The Jacobian contains non-finite",
                                         "values at iteration 1.\n")))
})

test_that("starting value 0.1", {
  expect_warning(
    expect_known_output(
      ret <- umf_solve_nl(0.1, f, jac, control = list(trace = TRUE)),
      "expected_output/singular2_report2.txt"),
    "NaNs produced"
  )
  expect_false(ret$solved)
  expect_true(is.nan(ret$fval))
  expect_true(ret$x != 0.1)
  expect_equal(ret$message, paste("Function value contains non-finite values",
            "(starting at index=1) at iteration 1\n"))
})
