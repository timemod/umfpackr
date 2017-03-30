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


xstart <- c(2,3)
ret <- umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                    control = list(trace = TRUE, silent = FALSE))

test_that("result of umf_solve is correct", {
    expect_true(ret$solved)
    expect_equal(ret$x, c(1,1))
    expect_true(sum(abs(ret$fval)) < 1e-8)
})
