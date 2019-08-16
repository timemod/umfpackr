library(umfpackr)
library(testthat)
library(methods)

rm(list = ls())

context("Test umf_solve_nl for a log function")

#x <- seq(0, 4, length.out = 100)[-1]
#plot(x, log(x), type = "l")
#lines(x,  0.5 * ( x - 1), type = "l", col = "red")
#lines(x,  2 * ( x - 1), type = "l", col = "green")

f <- function(x) {
  y <- numeric(2)
  y[1] <- log(x[1]) - 0.5 * (x[1] - 1)
  y[2] <- log(x[2]) - 2 * (x[2] - 1)
  return(y)
}

jacf <- function(x) {
  m <- matrix(0, 2, 2)
  m[1, 1] <- 1 / x[1] - 0.5
  m[2, 2] <- 1 / x[2] - 2
  return(as(m, "dgCMatrix"))
}

test_that("correct solution with appropriate starting point", {
  ret <- umf_solve_nl(c(1.2, 0.8), f, jacf,
                      control = list(silent = TRUE, trace = TRUE))
  expect_equal(ret$message, "ok")
  expect_equal(ret$x, c(1, 1))
})

test_that("linesearching", {

  # without linesearching, the function fails
  expect_warning({
    expect_output({
      ret <- umf_solve_nl(c(1.8, 1.1), f, jacf,
                     control = list(silent = TRUE, trace = TRUE))
    }, NA)
  }, "NaNs produced")
  expect_equal(ret$message,
           "Function value contains non-finite values (starting at index=1) at iteration 1\n")
  expect_false(ret$solved)

  # now use linesearching
  expect_silent(
    expect_warning(
      ret <- umf_solve_nl(c(1.8, 1.1), f, jacf,
                          control = list(silent = TRUE, trace = FALSE),
                          global = "cline")
    )
  )
  expect_equal(ret$message, "ok")
  expect_equal(ret$x, c(1, 1))
})

test_that("starting at singular point x = 2", {

  ret <- umf_solve_nl(c(2, 1.1),f, jacf, control = list(silent = TRUE))
  expect_equal(ret$message,
      paste("The Jacobian is (nearly) singular at iteration 1. The inverse",
            "condition is 0.\n"))
  expect_false(ret$solved)


  expect_silent(
    ret <- umf_solve_nl(c(2, 1.1), f, jacf, control = list(silent = TRUE,
                                                         trace = TRUE,
                                                 allow_singular = TRUE))
  )

  expect_equal(ret$message,
               paste("Relative step size smaller than xtol (1e-08)\n"))
  expect_false(ret$solved)
})





