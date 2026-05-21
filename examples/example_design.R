# End-to-end example: find all admissible two-stage designs
# for p0=0.10, p1=0.30, alpha=0.05, power=0.80, n_max=40.
#
# Run after installing the package:  Rscript examples/example_design.R

library(twostage)

cat("Searching for admissible designs...\n")
result <- admissible_designs(
  p0    = 0.10,
  p1    = 0.30,
  alpha = 0.05,
  power = 0.80,
  n_max = 40L
)
cat(sprintf("Found %d admissible designs on the 3D Pareto frontier.\n\n",
            nrow(result)))

display <- result[, c("n1", "n2", "n", "r1", "e1", "r",
                      "alpha_actual", "power_actual",
                      "en_null", "en_alt", "design_type")]
display$alpha_actual <- round(display$alpha_actual, 4)
display$power_actual <- round(display$power_actual, 4)
display$en_null      <- round(display$en_null, 2)
display$en_alt       <- round(display$en_alt, 2)
print(display, row.names = FALSE)

cat(sprintf("\nMinimax design: n = %d\n",
            result$n[grepl("minimax", result$design_type)]))
cat(sprintf("Optimal design: EN_null = %.2f\n",
            result$en_null[grepl("optimal", result$design_type)]))

cat("\n--- Simon two-stage designs (no interim efficacy stop) ---\n")
simon_result <- admissible_designs(
  p0    = 0.10,
  p1    = 0.30,
  alpha = 0.05,
  power = 0.80,
  n_max = 40L,
  simon = TRUE
)
cat(sprintf("Found %d admissible Simon designs.\n\n", nrow(simon_result)))

simon_display <- simon_result[, c("n1", "n2", "n", "r1", "e1", "r",
                                  "alpha_actual", "power_actual",
                                  "en_null", "en_alt", "design_type")]
simon_display$alpha_actual <- round(simon_display$alpha_actual, 4)
simon_display$power_actual <- round(simon_display$power_actual, 4)
simon_display$en_null      <- round(simon_display$en_null, 2)
simon_display$en_alt       <- round(simon_display$en_alt, 2)
print(simon_display, row.names = FALSE)
