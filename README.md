## umfpackr

Package `umfpackr` is an R interface to UMFPACK, which can be used
to solve linear equation using sparse matrix techniques. UMFPACK is part of the SuiteSparse package.
The R package `umfpackr` also implements a function to solve non-linear 
equations by a Newton method employing the sparse matrix method of UMFPACK.
The non-linear solvers included a simple cubic line search.

On Linux and other non-Windows system, the package employs library `libsuitesparse-dev`, which should
be installed on the system. For Ubuntu you can use `sudo apt-get install libsuitesparse-dev` to install this package.

On Windows, package `umfpackr` uses the source code of the SuiteSparse modules `UMFPACK` and `AMD`. This means that the METIS ordering
method is not possible.

## Documentation

[Reference manual](umfpackr.pdf)

[Vignette *UMFPACK Interface for R*](pkg/vignettes/UMFPACK_interface.pdf)

