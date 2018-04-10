# http://stackoverflow.com/questions/29134996/solving-nonlinear-equation-in-r

# wants to know if system has closed form solution
# I want to see how umfpackr behaves
#
# Based on test singular2.R of package nleqslv

set.seed(29)

library(testthat)
library(umfpackr)
library(Matrix)

context("Test option allow_singular for umfpack_solve_nl")

rm(list = ls())

f <- function(X, a, b, c1, c2, c3) {
  Y <- numeric(3)
  x <- X[1]
  y <- X[2]
  z <- X[3]
  Y[1] <- x + y - x * y - c1
  Y[2] <- x + z - x * z - c2
  Y[3] <- a * y + b * z - c3
  return(Y)
}

Jac <- function(X, a, b, c1, c2, c3) {
  J <- matrix(0,nrow=3,ncol=3)
  x <- X[1]
  y <- X[2]
  z <- X[3]

  J[1,1] <- 1-y
  J[2,1] <- 1-z
  J[3,1] <- 0
  J[1,2] <- 1-x
  J[2,2] <- 0
  J[3,2] <- a
  J[1,3] <- 0
  J[2,3] <- 1-x
  J[3,3] <- b

  return(as(J, "dgCMatrix"))
}

a <- 1
b <- 1
c1 <- 2
c2 <- 3
c3 <- 4

# exact solution
x <- (a * c1 + b * c2 - c3) / (a + b - c3)
y <- (b * c1 - b * c2 - c1 * c3 + c3) / (-a * c1 + a - b * c2 + b)
z <- (a * (c1 - c2) + (c2 - 1) * c3) / (a * (c1 - 1) + b * (c2-1))
xsol <- c(x, y ,z)

x_start <- c(1, 2 , 3)

test_that("without allow_singular", {
  msg <- "No convergence after 1 iterations"
  expect_output(z1 <- umf_solve_nl(x_start, f, Jac, a = a, b = b, c1 = c1,
                                   c2 = c2, c3 = c3,
                                   control = list(trace = FALSE)),
                msg)
  expect_false(z1$solved)
  expect_equal(z1$message,
            "The Jacobian is (nearly) singular. The inverse condition is 0.\n")
})

test_that("with allow_singular", {
  msg <- "Convergence after 3 iterations"
  expect_output(z1 <- umf_solve_nl(x_start, f, Jac, a = a, b = b, c1 = c1,
                       c2 = c2, c3 = c3,
                       control = list(trace = FALSE,
                                      allow_singular = TRUE)),
                msg)
  expect_true(z1$solved)
  expect_equal(z1$x, xsol)
})
