library(nleqslv)

dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}

jacdsln <- function(x) {
    n <- length(x)
    Df <- matrix(numeric(n*n),n,n)
    Df[1,1] <- 2*x[1]
    Df[1,2] <- 2*x[2]
    Df[2,1] <- exp(x[1]-1)
    Df[2,2] <- 3*x[2]^2

    Df
}


xstart <- c(2,3)
ret <- nleqslv(xstart, fn = dslnex, jac = jacdsln, method = "Newton", 
               control = list(trace = TRUE))
print(ret$x)
