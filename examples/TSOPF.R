library(Matrix)
library(umfpackr)

rm(list = ls())

# Example TSOPF_RS_b678_c2, downloaded from SuiteSparse Matrix Collection
# (https://sparse.tamu.edu).

# unzip data file. Note that the matrix file is zipped because
# GitHub does not accepts files with size > 100 MB
unzip("input/TSOPF_RS_b678_c2.mtx.zip")

mat_data <- read.table("TSOPF_RS_b678_c2.mtx"
                , comment.char = "%", row.names = NULL)

# remove the large file
file.remove("\TSOPF_RS_b678_c2.mtx")


nrow <- mat_data[1, 1]
ncol <- mat_data[1, 2]
nnz <- mat_data[1, 3]
mat_ata <- mat_data[-1, ]

sel <- mat_data[ , 3] != 0
mdat_data <- mat_data[sel, ]

rows <- mat_data[ , 1]
cols <- mat_data[ , 2]
values <- mat_data[ , 3]

m <- sparseMatrix(i = rows, j = cols, x = values, dims = c(nrow, ncol))


b_data <- read.table("input/TSOPF_RS_b678_c2_b.mtx"
                       , comment.char = "%", row.names = NULL)
row_b <- b_data[1, 1]
ncol_b <- b_data[1, 2]
nnz_b <- b_data[1, 3]
mat_ata <- b_data[-1, ]

sel <- b_data[ , 3] != 0
b_data <- b_data[sel, ]

rows_b <- b_data[ , 1]
cols_b <- b_data[ , 2]
values_b <- b_data[ , 3]

b <- sparseMatrix(i = rows_b, j = cols_b, x = values_b, dims = c(nrow, ncol))
b <- as.numeric(b[ , 1])

print(system.time(
  solve(m, b)
))

print(system.time(
  umf_solve(m, b)
))


