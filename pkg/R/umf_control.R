#' Argument `umf_control`
#'
#' Functions \code{\link{umf_solve}} and \code{\link{umf_solve_nl}} have an
#' argument `umf_control` that can be used to specify UMFPACK control parameters.
#' These control parameters are described in detail in the
#' \href{../doc/UMFPACK_UserGuide.pdf}{UMFPACK User Guide}. This guide describes
#' how the control options can be specified in C code. In the following paragraphs
#' is explained how the C code can be translated to an R expression passed to
#' argument `umf_control`.
#' \cr\cr
#' In C code, the control parameters are stored
#' in a numeric array `control` with elements specified with named numerical
#' constants. For example, element `UMFPACK_ORDERING` contains a control
#' parameter that specifies the ordering method.  The allowed values for this
#' parameters are specified with constants such as `UMFPACK_ORDERING_AMD`,
#' `UMFPACK_ORDERING_CHOLMOD` and `UMFPACK_ORDERING_METIS`.
#' For example, to use METIS, the specification in C code is:
#' ```
#' control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
#' ```
#' The corresponding R expression passed to argument `umf_control` is
#' ```
#' list(ORDERING = "METIS")
#' ```
#' Thus in the R expression the control names are without
#' prefix `UMFPACK_` and the possible values are specified with a character
#' giving the name of the constant in UMFPACK excluding the prefix
#' `UMFPACK_ORDERING_`.
#' \cr\cr
#' Some UMFPACK control parameters are specified with numerical values.
#' For example,
#' ```
#' list(SYM_PIVOT_TOLERANCE = 0.05)
#' ```
#' specifies the value of the symmetric pivot tolerance (the corresponding
#' C code would be
#' `control[UMFPACK_SYM_PIVOT_TOLERANCE] = 0.05`)
#' \cr\cr
#' Another example of an R expression for argument `umf_control`:
#' ```
#' list(STRATEGY = "UNSYMMETRIC", ORDERING = "METIS",
#'      SCALE = "NONE")
#' ```
#' @section Warning:
#'
#' On Windows, `umfpackr` does not use the complete SuiteSparse package
#' but onbly the `UMFPACK` and `AMD` modules. Therefore the only allowed
#' ordering option on Windows is `AMD`.
#' @seealso \code{\link{umf_solve}} and \code{\link{umf_solve_nl}}
#' @name umf_control
NULL
