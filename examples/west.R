library(umfpackr)
library(Matrix)


# Example west, downloaded from SuiteSparse Matrix Collection
# (https://sparse.tamu.edu).


x <- read.csv("input/west_acc.txt", header = FALSE)[ , 1]
i <- read.csv("input/west_icc.txt", header = FALSE)[ , 1] + 1
p <- read.csv("input/west_ccc.txt", header = FALSE)[ , 1]

n <- max(i)

mat <- sparseMatrix(i = i, p = p, x = x)

set.seed(1234545)
b <- runif(n)


library(microbenchmark)

t1 <- microbenchmark(
  x1 <- umf_solve(mat, b)
)
print(t1)

t2 <- microbenchmark(
  x2 <- solve(mat, b)
)
print(t2)

print(all.equal(x1, as.numeric(x2)))

t3 <- microbenchmark(
  x3 <- base::solve(as.matrix(mat), b)
)
print(t3)
