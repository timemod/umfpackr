#' @export
umf_solve_nl <- function(start, fn, jac, ..., control = list()) {

    control_ <- list(ftol = 1e-8, maxiter = 20, trace = TRUE,
                     cndtol = 1e-12, cnmtrx = 0.9)

    control_[names(control)] <- control

    fun <- function(x) {
        return(fn(x, ...))
    }
    jacob <- function(x) {
        return(jac(x, ...))
    }

    return(umf_solve_nl_(start, fun, jacob, control_))

}
