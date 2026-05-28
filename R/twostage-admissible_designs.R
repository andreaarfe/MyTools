#' Find all admissible two-stage designs
#'
#' Combines [twostage_find_feasible_designs()] and
#' [twostage_find_admissible_designs()] into a single call. Returns the
#' Pareto-optimal designs on three criteria —
#' maximum sample size (`n`), expected sample size under the null (`en_null`),
#' and expected sample size under the alternative (`en_alt`) — and labels four
#' designs of interest: the **minimax** (minimum `n`), the **optimal null**
#' (minimum `en_null`), the **optimal alt** (minimum `en_alt`), and the
#' **optimal mean** (minimum `0.5 * (en_null + en_alt)`) designs.
#'
#' @param p0          Numeric in (0, 1). Null response rate.
#' @param p1          Numeric in (0, 1). Alternative response rate (p1 > p0).
#' @param alpha       Numeric in (0, 1). Maximum allowable type I error.
#' @param power       Numeric in (0, 1). Minimum required power.
#' @param n_max       Integer. Maximum total sample size to consider (default 50).
#' @param n1_min      Integer. Minimum stage-1 sample size (default 1).
#' @param non_binding Logical. If `TRUE`, restrict to designs where `e1 >= r`,
#'   so that an interim efficacy declaration cannot be overturned at the final
#'   analysis. Use `TRUE` when the interim efficacy stopping rule is
#'   **non-binding** (the trial may continue to stage 2 even after `X1 >= e1`),
#'   so that the final decision remains coherent with the interim declaration.
#'   Default `FALSE`.
#' @param simon Logical. If `TRUE`, restrict to Simon two-stage designs,
#'   i.e. designs with no interim stopping for efficacy (`e1 = n1 + 1`).
#'   Default `FALSE`. When `simon = TRUE` the `non_binding` argument has no
#'   effect: because there is no interim efficacy declaration, the
#'   `e1 >= r` constraint is vacuously satisfied regardless of `r`.
#'
#' @return A `data.frame` with one row per admissible design and columns
#'   `n`, `n1`, `n2`, `r1`, `e1`, `r`, `alpha_actual`,
#'   `power_actual`, `en_null`, `en_alt`, `is_non_binding`, and `design_type`.
#'   The `design_type` column is a character vector with values:
#'   \describe{
#'     \item{`"minimax"`}{The admissible design with the uniquely smallest `n`.
#'       When several designs share the minimum `n`, the one with the smallest
#'       `en_null` (then `en_alt`, then `en_mean`) is chosen, and every secondary
#'       criterion it also wins within that tied set is appended to the label
#'       (e.g. `"minimax, min en_null"`, `"minimax, min en_null, min en_alt, min en_mean"`).}
#'     \item{`"optimal null"`}{The admissible design with the smallest `en_null`
#'       (ties broken by `n`).}
#'     \item{`"optimal alt"`}{The admissible design with the smallest `en_alt`
#'       (ties broken by `n`).}
#'     \item{`"optimal mean"`}{The admissible design with the smallest
#'       `0.5 * (en_null + en_alt)` (ties broken by `n`).}
#'     \item{`""`}{All other admissible designs.}
#'   }
#'   A design satisfying several criteria carries a comma-joined label
#'   (e.g. `"minimax, optimal null"`).
#'   Returns an empty `data.frame` when no feasible designs are found.
#'
#' @examples
#' twostage_admissible_designs(p0 = 0.10, p1 = 0.30, alpha = 0.05, power = 0.80,
#'                             n_max = 40L)
#' @export
twostage_admissible_designs <- function(p0, p1, alpha, power,
                                        n_max = 50L, n1_min = 1L,
                                        non_binding = FALSE,
                                        simon = FALSE) {
  feasible  <- twostage_find_feasible_designs(p0, p1, alpha, power,
                                              n_max, n1_min, non_binding, simon)
  if (nrow(feasible) == 0L) return(feasible)

  admissible <- twostage_find_admissible_designs(feasible)
  if (nrow(admissible) == 0L) return(admissible)

  # minimax: smallest n; when tied, narrow sequentially by en_null, en_alt,
  # en_mean; label reflects every secondary criterion the chosen design wins
  # within the full minimax candidate set.
  all_mm   <- which(admissible$n == min(admissible$n))
  mm_idxs  <- all_mm
  mm_idxs  <- mm_idxs[admissible$en_null[mm_idxs] == min(admissible$en_null[mm_idxs])]
  mm_idxs  <- mm_idxs[admissible$en_alt [mm_idxs] == min(admissible$en_alt [mm_idxs])]
  en_mean_mm <- 0.5 * (admissible$en_null[mm_idxs] + admissible$en_alt[mm_idxs])
  mm_idxs  <- mm_idxs[en_mean_mm == min(en_mean_mm)]
  mm_idx   <- mm_idxs[[1L]]
  mm_parts <- "minimax"
  if (length(all_mm) > 1L) {
    en_mean_all <- 0.5 * (admissible$en_null[all_mm] + admissible$en_alt[all_mm])
    if (admissible$en_null[mm_idx] == min(admissible$en_null[all_mm]))
      mm_parts <- c(mm_parts, "min en_null")
    if (admissible$en_alt[mm_idx]  == min(admissible$en_alt[all_mm]))
      mm_parts <- c(mm_parts, "min en_alt")
    if (0.5 * (admissible$en_null[mm_idx] + admissible$en_alt[mm_idx]) == min(en_mean_all))
      mm_parts <- c(mm_parts, "min en_mean")
  }
  mm_label <- paste(mm_parts, collapse = ", ")

  # optimal null: smallest en_null, tie-break by n
  on_idx <- which(admissible$en_null == min(admissible$en_null))
  on_idx <- on_idx[which.min(admissible$n[on_idx])]

  # optimal alt: smallest en_alt, tie-break by n
  oa_idx <- which(admissible$en_alt == min(admissible$en_alt))
  oa_idx <- oa_idx[which.min(admissible$n[oa_idx])]

  # optimal mean: smallest 0.5 * (en_null + en_alt), tie-break by n
  en_mean <- 0.5 * (admissible$en_null + admissible$en_alt)
  om_idx  <- which(en_mean == min(en_mean))
  om_idx  <- om_idx[which.min(admissible$n[om_idx])]

  labels <- vector("list", nrow(admissible))
  labels[[mm_idx]] <- c(labels[[mm_idx]], mm_label)
  labels[[on_idx]] <- c(labels[[on_idx]], "optimal null")
  labels[[oa_idx]] <- c(labels[[oa_idx]], "optimal alt")
  labels[[om_idx]] <- c(labels[[om_idx]], "optimal mean")

  admissible$design_type <- vapply(
    labels,
    function(x) if (is.null(x)) "" else paste(x, collapse = ", "),
    character(1))
  admissible
}
