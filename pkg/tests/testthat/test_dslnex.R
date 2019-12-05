library(umfpackr)
library(testthat)
library(methods)

context("Test umf_solve_nl with dslnex mnodel")

dslnex <- function(x, c) {
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - c
  y[2] <- exp(x[1]-1) + x[2]^3 - c
  y
}

jacdsln <- function(x, c) {
  n <- length(x)
  Df <- matrix(numeric(n*n),n,n)
  Df[1,1] <- 2*x[1]
  Df[1,2] <- 2*x[2]
  Df[2,1] <- exp(x[1]-1)
  Df[2,2] <- 3*x[2]^2

  return(as(Df, "dgCMatrix"))
}

test_that("result of umf_solve is correct", {
  xstart <- c(2, 3)
  ret <- umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                      control = list(trace = FALSE, silent = TRUE))
  expect_true(ret$solved)
  expect_equal(ret$x, c(1, 1))
  expect_true(sum(abs(ret$fval)) < 1e-8)

  expect_error(umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                            umf_control = list(SCALE = "NONE")),
               "When row scaling is applied, then umfpack control SCALE should not be \"SCALE_NONE\".")
})

test_that("result of umf_solve with cublic line search is correct", {
  xstart <- c(2, 0.5)

  ret <- umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                      control = list(trace = FALSE, silent = TRUE),
                      global = "cline", scaling = "row",
                      umf_control = list(SCALE = "SUM"))
  expect_true(ret$solved)
  expect_equal(ret$x, c(1, 1))
  expect_true(sum(abs(ret$fval)) < 1e-8)

  expect_warning(umf_solve_nl(xstart, dslnex, jacdsln, c = 2, scaling = "no",
                              umf_control = list(SCALE = "SUM")),
                 "Umpack control SCALE is ignored if row scaling is disabled.")
})

test_that("result of umf_solve without cublic line search is correct", {
  xstart <- c(2, 0.5)
  if (.Platform$OS.type == "windows") {
    expect_error(umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                              umf_control = list(ORDERING = "METIS")),
                 "On Windows, only the AMD ordering method possible.")
    ret <- umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                        control = list(trace = FALSE, silent = TRUE),
                        global = "no", umf_control = list(ORDERING = "AMD"))
  } else {
    ret <- umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                        control = list(trace = FALSE, silent = TRUE),
                        global = "no", umf_control = list(ORDERING = "METIS"))
  }
  expect_false(ret$solved)
})



