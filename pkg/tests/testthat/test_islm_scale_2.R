library(umfpackr)
library(testthat)
library(nleqslv)
library(stringr)

rm(list = ls())

context("islm model (row scaling)")

source("tools/read_rconds.R")

#
# parameter values
#

c0 <- 100
c1 <- 0.32
c2 <- -20
c3 <- 1

i0 <- 100
i1 <- 0.08
i2 <- -40
i3 <- 1.5

m0 <- 75
m1 <- 0.23
m2 <- -35
m3 <- 1.5

t0 <- -15
t1 <- 0.22



fun <- function(x) {
  x_list <- as.list(x)
  list2env(as.list(x), environment())

  res <- numeric(7)
  res[1] <- (y  - (c + i + g)) * scale
  res[2] <- yd - (y - t)
  res[3] <- t - (t0  + t1 * y)
  res[4] <- c - (c0  + c1 * yd + c2 * r + c3 * r^2)
  res[5] <- i - (i0  + i1 * y + i2 * r + i3 * r^2)
  res[6] <- md - (m0 + m1 * y + m2 * r + m3 * r^2)
  res[7] <- md - ms
  return(res)
}

jac <- function(x) {
  x_list <- as.list(x)
  list2env(as.list(x), environment())

  n <- 7
  j <- matrix(numeric(n * n), n ,n)
  colnames(j) <- names(x)

  j[1, "y"] <- 1 * scale
  j[1, "c"] <- -1 * scale
  j[1, "i"] <- -1 * scale

  j[2, "yd"] <- 1
  j[2, "y"] <- -1
  j[2, "t"] <- 1

  j[3, "t"] <- 1
  j[3, "y"] <- -t1

  j[4, "c"] <- 1
  j[4, "yd"] <- -c1
  j[4, "r"] <- - (c2 + 2 * c3 * r)

  j[5, "i"] <- 1
  j[5, "y"] <- -i1
  j[5, "r"] <- - (i2 + 2 * i3 * r)

  j[6, "md"] <- 1
  j[6, "y"] <- -m1
  j[6, "r"] <- - (m2 + 2 * m3 * r)

  j[7, "md"] <- 1

  return(as(j, "dgCMatrix"))
}



g <- 240
ms <- 100
xstart<- c(y = 460, yd = 371, t = 90, c = 170, i = 40, md = 100, r = 2.5)

scale <- 1
result_scale1  <- umf_solve_nl(xstart, fun, jac,
                               control = list(trace = FALSE, silent = TRUE))


#
# scale 1e32
#
scale <- 1e32

test_that("scale 1e32)", {

  report1 <- capture.output(
    result1 <- umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE)))
  expect_true(result1$solved)
  expect_equal(result1$x, result_scale1$x)
  rconds <- read_rconds(report1)
  expect_true(all(rconds > 1e-3))

  report1a_error <- capture_output(
    result1a_error <- umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE),
                             scaling = "no"))
  msg <- paste("The inverse condition of the jacobian is smaller",
               "than cnd_tol")
  expect_true(grepl(msg, report1a_error))
  expect_true(grepl(paste0("^", msg), result1a_error$message))

  #
  # test cnd_method == "condest"
  #
  set.seed(123) # specify random seed, needed for the condest method
  report1b <- capture.output(
    result1b <- umf_solve_nl(xstart, fun, jac,
                             control = list(trace = TRUE,
                                            cnd_method = "condest")))
  expect_true(result1b$solved)
  expect_equal(result1b$x, result_scale1$x)
  rconds <- read_rconds(report1b)
  expect_true(all(rconds > 1e-3))

  cond_iter1 <- as.numeric(sub("^\\s*\\d+\\s*(\\S+).+", "\\1",
                              report1b, perl = TRUE)[6])
  set.seed(123)
  j <- jac(xstart)
  j <- j / Matrix::rowSums(abs(j))
  expected <- 1 / Matrix::condest(j)$est
  expect_equal(cond_iter1, expected, tolerance = 1e-3)

  #
  # test cnd_method == "kappa"
  #
  report1c <- capture.output(
    result1c <- umf_solve_nl(xstart, fun, jac,
                             control = list(trace = TRUE,
                                            cnd_method = "kappa")))
  expect_true(result1c$solved)
  expect_equal(result1c$x, result_scale1$x)
  rconds <- read_rconds(report1c)
  expect_true(all(rconds > 1e-3))
  cond_iter1 <- as.numeric(sub("^\\s*\\d+\\s*(\\S+).+", "\\1",
                               report1c, perl = TRUE)[6])
  j <- jac(xstart)
  j <- j / Matrix::rowSums(abs(j))
  expected <- 1 / kappa(j, exact = TRUE)
  expect_true(abs(cond_iter1 - expected) < 1e-4)

  #
  # test column  scaling
  #
  report2_error <- capture_output(
    result2_error <- umf_solve_nl(xstart, fun, jac,
                                  control = list(trace = TRUE),
                                  scaling = "col")
  )
  expect_false(result2_error$solved)
  msg <- paste("The inverse condition of the jacobian is smaller",
               "than cnd_tol")
  expect_true(grepl(msg, report2_error))
  expect_true(grepl(paste0("^", msg), result2_error$message))


  report2 <- capture.output(
    result2 <- umf_solve_nl(xstart, fun, jac,
                           control = list(trace = TRUE,
                                          cnd_method = "kappa",
                                          cnd_tol = 0), scaling = "col"))
  expect_true(result2$solved)
  expect_equal(result2$x, result_scale1$x)
  rconds <- read_rconds(report2)
  expect_true(all(rconds < 1e-32))

  # no scaling:
  report3 <- capture.output(
    result3 <- umf_solve_nl(xstart, fun, jac,
                           control = list(trace = TRUE,
                                          cnd_tol = -999), scaling = "no"))
  expect_true(result3$solved)
  expect_equal(result3$x, result_scale1$x)
  rconds <- read_rconds(report3)
  expect_true(all(rconds < 1e-32))

  # allow_singular does not help here
  report4_error <- capture_output(
    result4_error <- umf_solve_nl(xstart, fun, jac,
                            control = list(trace = TRUE,
                                           allow_singular = TRUE),
                            scaling = "no"))
  expect_false(result4_error$solved)
  msg <- "Relative step size smaller than"
  expect_true(grepl(msg, report4_error))
  expect_true(grepl(paste0("^", msg), result4_error$message))
})
