#' Filter designs to the Pareto frontier on sample size and expected sample sizes
#'
#' Retains only admissible (Pareto-optimal) designs with respect to three
#' criteria: maximum sample size (`n`), expected sample size under the null
#' (`en_null`), and expected sample size under the alternative (`en_alt`).
#' A design is admissible if no other design is weakly better on all three
#' axes with at least one strict improvement.
#'
#' @param designs A `data.frame` as returned by [find_feasible_designs()],
#'   containing at least columns `n`, `en_null`, and `en_alt`.
#'
#' @return A `data.frame` of admissible designs sorted by `en_null` ascending.
#'   Returns `designs` unchanged (empty) when the input has zero rows.
#' @export
find_admissible_designs <- function(designs) {
  if (nrow(designs) == 0L) return(designs)

  nm  <- designs$n
  en0 <- designs$en_null
  en1 <- designs$en_alt
  m   <- nrow(designs)
  keep <- rep(TRUE, m)

  for (i in seq_len(m)) {
    if (!keep[i]) next
    # Check whether any j dominates i on all three axes
    dominated <- any(
      keep & seq_len(m) != i &
        nm  <= nm[i]  & en0 <= en0[i] & en1 <= en1[i] &
        (nm  < nm[i]  | en0 < en0[i]  | en1 < en1[i])
    )
    if (dominated) keep[i] <- FALSE
  }

  designs[keep, ][order(designs$en_null[keep]), ]
}
