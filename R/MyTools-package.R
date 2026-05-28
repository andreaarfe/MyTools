#' @description
#' A collection of R utilities for statistical and clinical-trial methodology.
#' The current release provides the **twostage** module (functions prefixed
#' `twostage_`), which enumerates admissible two-stage single-arm clinical
#' trial designs for a binary endpoint. The [twostage_admissible_designs()]
#' entry point returns the 3D Pareto frontier on `(n, en_null, en_alt)`; the
#' `non_binding = TRUE` option enforces the constraint `e1 >= r` so that an
#' interim efficacy declaration cannot be overturned by the final analysis,
#' and `simon = TRUE` restricts the search to Simon two-stage designs (no
#' interim efficacy stop).
#' @keywords internal
#' @seealso \code{vignette("twostage", package = "MyTools")}
"_PACKAGE"

#' @importFrom Rcpp sourceCpp
#' @useDynLib MyTools, .registration = TRUE
NULL
