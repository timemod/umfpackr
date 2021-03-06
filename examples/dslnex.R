library(umfpackr)

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

#xstart <- c(2,0.5)
xstart <- c(2, 0.5)

print(dslnex(xstart, c = 2))

x <- umf_solve_nl(xstart, dslnex, jacdsln, c = 2,
                  control = list(trace = TRUE, maxiter = 10, xtol = 1e-8),
                  global = "cline")

#ret <- nleqslv(xstart, fn = dslnex, jac = jacdsln, method = "Newton",
#               control = list(trace = TRUE))
#print(ret$x)
