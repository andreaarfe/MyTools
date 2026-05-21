#' @description
#' Enumerates all admissible two-stage single-arm clinical trial designs for a
#' binary endpoint. Supports futility and efficacy stopping at an interim
#' analysis. The `irrevocable = TRUE` option (the default) enforces the
#' constraint `e1 >= r`, ensuring that an interim efficacy declaration cannot
#' be overturned by the final analysis. This should be used when the interim
#' efficacy stopping rule is **non-binding** — i.e. the trial may continue to
#' stage 2 even after declaring interim success — so that the final decision
#' remains coherent with the interim declaration.
#' @keywords internal
#' @seealso \code{vignette("twostage", package = "twostage")}
"_PACKAGE"

#' @importFrom Rcpp sourceCpp
#' @useDynLib twostage, .registration = TRUE
NULL
