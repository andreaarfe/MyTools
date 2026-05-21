library(testthat)

# Resolve the package root regardless of working directory.
.find_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) == 1L) {
    return(normalizePath(file.path(dirname(file_arg), "..")))
  }
  # Fallback for testthat::test_file() and similar callers
  if (file.exists("R/distributions.R")) return(normalizePath("."))
  if (file.exists("../R/distributions.R")) return(normalizePath(".."))
  stop("Cannot locate package root")
}
.root <- .find_root()

source(file.path(.root, "R", "distributions.R"))
source(file.path(.root, "R", "design.R"))
source(file.path(.root, "R", "search.R"))
source(file.path(.root, "R", "admissibility.R"))

# ---------------------------------------------------------------------------
# distributions.R
# ---------------------------------------------------------------------------
test_that("binom_pmf matches dbinom", {
  expect_equal(binom_pmf(3, 10, 0.3), dbinom(3, 10, 0.3))
})

test_that("binom_cdf matches pbinom", {
  expect_equal(binom_cdf(4, 10, 0.2), pbinom(4, 10, 0.2))
})

test_that("binom_upper = P(X >= k)", {
  expect_equal(binom_upper(5, 10, 0.4), 1 - pbinom(4, 10, 0.4))
})

# ---------------------------------------------------------------------------
# design.R — evaluate_design
# ---------------------------------------------------------------------------

# With r1 = -1 (no futility) and e1 = n1+1 (no interim efficacy), the design
# collapses to a one-stage test: success iff X_total >= r.
# P(success) = P(Binom(n, p) >= r).
test_that("one-stage equivalence when interim stops disabled", {
  n1 <- 10L; n <- 20L; r <- 6L
  d  <- evaluate_design(n1, n, r1 = -1L, e1 = n1 + 1L, r = r, p0 = 0.2, p1 = 0.5)
  expect_equal(d$alpha_actual, pbinom(r - 1L, n, 0.2, lower.tail = FALSE),
               tolerance = 1e-10)
  expect_equal(d$power_actual, pbinom(r - 1L, n, 0.5, lower.tail = FALSE),
               tolerance = 1e-10)
  # EN must equal n when trial never stops early
  expect_equal(d$en_null, n, tolerance = 1e-10)
  expect_equal(d$en_alt,  n, tolerance = 1e-10)
})

test_that("is_irrevocable is TRUE iff e1 >= r", {
  for (r in 3:7) {
    for (e1 in 1:8) {
      d <- evaluate_design(n1 = 8L, n = 15L, r1 = -1L, e1 = e1, r = r,
                           p0 = 0.2, p1 = 0.5)
      expect_equal(d$is_irrevocable, e1 >= r,
                   label = sprintf("e1=%d, r=%d", e1, r))
    }
  }
})

test_that("EN is bounded above by n", {
  d <- evaluate_design(n1 = 10L, n = 20L, r1 = 5L, e1 = 8L, r = 8L,
                       p0 = 0.2, p1 = 0.5)
  expect_lte(d$en_null, 20)
  expect_lte(d$en_alt,  20)
})

test_that("success probability is 0 when r exceeds total sample size", {
  d <- evaluate_design(n1 = 5L, n = 10L, r1 = -1L, e1 = 6L, r = 11L,
                       p0 = 0.2, p1 = 0.5)
  expect_equal(d$alpha_actual, 0, tolerance = 1e-15)
  expect_equal(d$power_actual, 0, tolerance = 1e-15)
})

test_that("empty continuation region (r1 = e1 - 1) gives no stage-2 contribution", {
  # When r1 = e1 - 1, every outcome triggers either futility or efficacy at
  # interim, so the trial always stops at stage 1: P(success) = P(X1 >= e1)
  # and EN = n1.
  n1 <- 10L; n <- 20L; e1 <- 6L; r1 <- e1 - 1L; r <- 6L
  for (p in c(0.1, 0.3, 0.5, 0.8)) {
    d <- evaluate_design(n1, n, r1, e1, r, p0 = p, p1 = p)
    p_eff <- pbinom(e1 - 1L, n1, p, lower.tail = FALSE)
    expect_equal(d$alpha_actual, p_eff, tolerance = 1e-10,
                 label = sprintf("p=%.2f, alpha", p))
    expect_equal(d$power_actual, p_eff, tolerance = 1e-10,
                 label = sprintf("p=%.2f, power", p))
    expect_equal(d$en_null, n1, tolerance = 1e-10,
                 label = sprintf("p=%.2f, EN", p))
  }
})

test_that("success probability matches manual computation in continuation region", {
  # Hand-checked example.
  n1 <- 5L; n <- 10L; r1 <- 1L; e1 <- 4L; r <- 5L; p <- 0.3
  d <- evaluate_design(n1, n, r1, e1, r, p0 = p, p1 = p)

  # P(efficacy at stage 1) = P(X1 >= 4)
  p_eff1 <- pbinom(3L, n1, p, lower.tail = FALSE)
  # Continuation region: X1 in {2, 3}; need X2 >= r - x1
  p_cont <- dbinom(2L, n1, p) * pbinom(2L, 5L, p, lower.tail = FALSE) +
            dbinom(3L, n1, p) * pbinom(1L, 5L, p, lower.tail = FALSE)
  expected <- p_eff1 + p_cont
  expect_equal(d$alpha_actual, expected, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# search.R — find_feasible_designs
# ---------------------------------------------------------------------------

test_that("all returned designs satisfy constraints", {
  alpha <- 0.05; pw <- 0.80
  feasible <- find_feasible_designs(p0 = 0.1, p1 = 0.3, alpha = alpha,
                                    power = pw, n_max = 25L)
  expect_gt(nrow(feasible), 0)
  expect_true(all(feasible$alpha_actual <= alpha + 1e-9))
  expect_true(all(feasible$power_actual >= pw     - 1e-9))
  expect_true(all(feasible$is_irrevocable))
  expect_true(all(feasible$r1 < feasible$e1))
})

# ---------------------------------------------------------------------------
# admissibility.R — find_admissible_designs
# ---------------------------------------------------------------------------

test_that("admissible set is a subset of feasible designs", {
  feasible   <- find_feasible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                                      power = 0.80, n_max = 25L)
  admissible <- find_admissible_designs(feasible)
  expect_lte(nrow(admissible), nrow(feasible))
  keys_adm <- paste(admissible$n1, admissible$n, admissible$r1,
                    admissible$e1, admissible$r)
  keys_fea <- paste(feasible$n1,   feasible$n,   feasible$r1,
                    feasible$e1,   feasible$r)
  expect_true(all(keys_adm %in% keys_fea))
})

test_that("no admissible design is dominated by another admissible design", {
  feasible   <- find_feasible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                                      power = 0.80, n_max = 25L)
  admissible <- find_admissible_designs(feasible)
  m   <- nrow(admissible)
  en0 <- admissible$en_null
  en1 <- admissible$en_alt
  for (i in seq_len(m)) {
    dominated <- any(
      seq_len(m) != i &
        en0 <= en0[i] & en1 <= en1[i] &
        (en0 < en0[i] | en1 < en1[i])
    )
    expect_false(dominated, label = sprintf("design %d is dominated", i))
  }
})

test_that("admissible set from synthetic designs has correct Pareto frontier", {
  # D1 dominated by D2 on both axes; D3 on frontier (trades EN0 for EN1).
  make_row <- function(en0, en1) {
    data.frame(n1=10, n=20, n2=10, r1=2, e1=8, r=8,
               p0=0.1, p1=0.3,
               alpha_actual=0.04, power_actual=0.82,
               en_null=en0, en_alt=en1, is_irrevocable=TRUE)
  }
  designs <- rbind(make_row(15, 12),   # D1 — dominated by D2
                   make_row(14, 11),   # D2 — dominates D1
                   make_row(13, 13))   # D3 — on frontier
  adm <- find_admissible_designs(designs)
  expect_equal(nrow(adm), 2L)
  expect_true(all(adm$en_null %in% c(14, 13)))
})
