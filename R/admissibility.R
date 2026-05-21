# Filter a data.frame of feasible designs to the Pareto frontier on
# (en_null, en_alt): a design is admissible if no other feasible design has
# EN_null <= its EN_null AND EN_alt <= its EN_alt, with at least one strict.
#
# Returns a data.frame of admissible designs sorted by en_null ascending.
find_admissible_designs <- function(designs) {
  if (nrow(designs) == 0L) return(designs)

  en0 <- designs$en_null
  en1 <- designs$en_alt
  m   <- nrow(designs)
  keep <- rep(TRUE, m)

  for (i in seq_len(m)) {
    if (!keep[i]) next
    # Check whether any j dominates i
    dominated <- any(
      keep & seq_len(m) != i &
        en0 <= en0[i] & en1 <= en1[i] &
        (en0 < en0[i] | en1 < en1[i])
    )
    if (dominated) keep[i] <- FALSE
  }

  designs[keep, ][order(designs$en_null[keep]), ]
}
