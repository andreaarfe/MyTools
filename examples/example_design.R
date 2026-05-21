# End-to-end example: find all admissible two-stage designs
# for p0=0.10, p1=0.30, alpha=0.05, power=0.80, n_max=40.
#
# Run from the package root:  Rscript examples/example_design.R

source("R/distributions.R")
source("R/design.R")
source("R/search.R")
source("R/admissibility.R")

cat("Searching for feasible designs...\n")
feasible <- find_feasible_designs(
  p0    = 0.10,
  p1    = 0.30,
  alpha = 0.05,
  power = 0.80,
  n_max = 40L
)
cat(sprintf("Found %d feasible designs.\n\n", nrow(feasible)))

admissible <- find_admissible_designs(feasible)
cat(sprintf("Admissible designs (%d on Pareto frontier of EN_null vs EN_alt):\n",
            nrow(admissible)))

display <- admissible[, c("n1", "n2", "n", "r1", "e1", "r",
                          "alpha_actual", "power_actual",
                          "en_null", "en_alt")]
display$alpha_actual <- round(display$alpha_actual, 4)
display$power_actual <- round(display$power_actual, 4)
display$en_null      <- round(display$en_null, 2)
display$en_alt       <- round(display$en_alt, 2)
print(display, row.names = FALSE)

# Sanity check: all returned designs must satisfy irrevocability
stopifnot(all(admissible$is_irrevocable))
cat("\nAll admissible designs satisfy the irrevocability constraint (e1 >= r).\n")
