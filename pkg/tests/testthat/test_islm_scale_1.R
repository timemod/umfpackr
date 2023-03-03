library(umfpackr)
library(testthat)
library(nleqslv)
library(stringr)

rm(list = ls())


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
  res[1] <- y  - (c + i + g)
  res[2] <- yd - (y - t)
  res[3] <- t - (t0 * scale + t1 * y)
  res[4] <- c - (c0 * scale + c1 * yd + (c2 * r + c3 * r^2) * scale)
  res[5] <- i - (i0 * scale + i1 * y + (i2 * r + i3 * r^2) * scale)
  res[6] <- md - (m0 * scale + m1 * y + (m2 * r + m3 * r^2) * scale)
  res[7] <- md - ms
  return(res / scale)
}

jac <- function(x) {
  x_list <- as.list(x)
  list2env(as.list(x), environment())

  n <- 7
  j <- matrix(numeric(n * n), n ,n)
  colnames(j) <- names(x)

  j[1, "y"] <- 1
  j[1, "c"] <- -1
  j[1, "i"] <- -1

  j[2, "yd"] <- 1
  j[2, "y"] <- -1
  j[2, "t"] <- 1

  j[3, "t"] <- 1
  j[3, "y"] <- -t1

  j[4, "c"] <- 1
  j[4, "yd"] <- -c1
  j[4, "r"] <- - (c2 + 2 * c3 * r) * scale

  j[5, "i"] <- 1
  j[5, "y"] <- -i1
  j[5, "r"] <- - (i2 + 2 * i3 * r) * scale

  j[6, "md"] <- 1
  j[6, "y"] <- -m1
  j[6, "r"] <- - (m2 + 2 * m3 * r) * scale

  j[7, "md"] <- 1

  j <- as(j / scale, "dgCMatrix")

  return(j)
}



g_scale_1 <- 240
ms_scale_1 <- 100
xstart_scale_1 <- c(y = 460, yd = 371, t = 90, c = 170, i = 40, md = 100,
                    r = 2.5)

#
# scale 1
#
scale <- 1
g <- g_scale_1
ms <- ms_scale_1
xstart <- xstart_scale_1

result_scale1  <- umf_solve_nl(xstart, fun, jac,
                        control = list(trace = FALSE, silent = TRUE))

test_that("check results scale 1", {

  # check results
  expect_true(result_scale1$solved)
  expect_known_value(result_scale1$x, file = "expected_output/islm_scale_result1.rds")


  # compare result with nleqslv result
  jac_nleqslv <-  function(x) {
    j <- jac(x)
    return(as.matrix(j))
  }

  result_nleqslv <- nleqslv(xstart, fun, jac_nleqslv,
                            control = list(trace = FALSE, ftol = 1e-8),
                            method = "Newton")

  expect_equal(result_scale1$x, result_nleqslv$x)
})


#
# scale 1e32
#
scale <- 1e32
g <- g_scale_1 * scale
ms <- ms_scale_1 * scale
xstart <- xstart_scale_1
xstart[1:6] <- xstart[1:6] * scale

test_that("scale 1e-12)", {

  # first test with default value of control parameter cnd_tol
  report_error <- capture_output(
    result_error <- umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE))
  )
  expect_false(result_error$solved)
  msg <- paste("The inverse condition of the jacobian is smaller",
               "than cnd_tol")
  expect_true(grepl(msg, report_error))
  expect_true(grepl(paste0("^", msg), result_error$message))

  report1 <- capture.output(
    result1 <- umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE,
                                                             cnd_tol = 0))
  )

  expect_true(result1$solved)

  expected_result <- result_scale1$x
  expected_result[1:6] <- scale * expected_result[1:6]
  expect_equal(result1$x, expected_result)

  rconds <- read_rconds(report1)
  expect_true(all(rconds < 1e-32))

  #
  # column scaling
  #
  report2 <- capture.output(
    result2 <- umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE),
                            scaling = "col"))
  expect_true(result2$solved)

  expected_result <- result_scale1$x
  expected_result[1:6] <- scale * expected_result[1:6]
  expect_equal(result2$x, expected_result)

  rconds <- read_rconds(report2)
  expect_true(all(rconds > 1e-3))


  # now with cnd_method == "kappa":
  report2a <- capture.output(
    result2a <- umf_solve_nl(xstart, fun, jac,
                            control = list(trace = TRUE, cnd_method = "kappa"),
                            scaling = "col"))
  expect_true(result2a$solved)

  expected_result <- result_scale1$x
  expected_result[1:6] <- scale * expected_result[1:6]
  expect_equal(result2a$x, expected_result)
  rconds <- read_rconds(report2a)
  expect_true(all(rconds > 1e-3))
  cond_iter1 <- as.numeric(sub("^\\s*\\d+\\s*(\\S+).+", "\\1",
                               report2a, perl = TRUE)[6])
  j <- jac(xstart)
  j <- j / rep(Matrix::colSums(abs(j)), each = nrow(j))
  expected <- 1 / kappa(j, exact = TRUE)
  expect_true(abs(cond_iter1 - expected) < 1e-3)

  report3 <- capture.output(
    result3 <- umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE,
                                                             cnd_tol = -1),
                           scaling = "no"))
  expect_true(result3$solved)

  expected_result <- result_scale1$x
  expected_result[1:6] <- scale * expected_result[1:6]
  expect_equal(result3$x, expected_result)

  rconds <- read_rconds(report3)
  expect_true(all(rconds < 1e-32))
})


test_that("errors", {
  msg <- "Control option 'cnd_tol' should be a numeric of length 1"
  expect_error(umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE,
                                                cnd_tol = "xxx")),
                            msg)
  expect_error(umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE,
                                                             cnd_tol = c(1,2))),
               msg)
  msg <- "Control option 'cnd_tol' should be a numeric of length 1"
  expect_error(umf_solve_nl(xstart, fun, jac, control = list(trace = TRUE,
                                                             cnd_tol = "xxx")),
               msg)

  msg <-  paste("Allowed values for control option 'cnd_method' are 'umfpack',",
                "'condest' and 'kappa'\\.")
  expect_error(umf_solve_nl(xstart, fun, jac, control = list(cnd_method = "xxx")),
               msg)

  msg <- "Control option 'cnd_method' should be a character of length 1"
  expect_error(umf_solve_nl(xstart, fun, jac,
                            control = list(cnd_method = c("kappa", "x"))),
               msg)
  expect_error(umf_solve_nl(xstart, fun, jac,
                            control = list(cnd_method = 2)),
               msg)

  msg <- "Unknown control options xxx."
  expect_error(umf_solve_nl(xstart, fun, jac,
                            control = list(xxx = 2)),
               msg)

})
