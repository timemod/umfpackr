library(umfpackr)

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

g <- 240
ms <- 100

fun <- function(x) {
  x_list <- as.list(x)
  list2env(as.list(x), environment())

  res <- numeric(7)
  res[1] <- y  - (c + i + g)
  res[2] <- yd - (y - t)
  res[3] <- t - (t0 + t1 * y)
  res[4] <- c - (c0 + c1 * yd + c2 * r + c3 * r^2)
  res[5] <- i - (i0 + i1 * y + i2 * r + i3 * r^2)
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


xstart <- c(y = 460, yd = 371, t = 90, c = 170, i = 40, md = 100,
            r = 2.5)

ret <- umf_solve_nl(xstart, fun, jac,
                    control = list(trace = TRUE, silent = FALSE))

print(ret$x)
