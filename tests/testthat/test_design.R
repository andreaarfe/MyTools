library(testthat)

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

test_that("all returned designs satisfy constraints (irrevocable = TRUE)", {
  alpha <- 0.05; pw <- 0.80
  feasible <- find_feasible_designs(p0 = 0.1, p1 = 0.3, alpha = alpha,
                                    power = pw, n_max = 25L,
                                    irrevocable = TRUE)
  expect_gt(nrow(feasible), 0)
  expect_true(all(feasible$alpha_actual <= alpha + 1e-9))
  expect_true(all(feasible$power_actual >= pw     - 1e-9))
  expect_true(all(feasible$is_irrevocable))
  expect_true(all(feasible$r1 < feasible$e1))
})

test_that("irrevocable=FALSE returns designs with e1 < r", {
  feasible <- find_feasible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                                    power = 0.80, n_max = 25L,
                                    irrevocable = FALSE)
  expect_gt(nrow(feasible), 0)
  # Must include at least some non-irrevocable designs
  expect_true(any(!feasible$is_irrevocable))
  # Error and power constraints still hold
  expect_true(all(feasible$alpha_actual <= 0.05 + 1e-9))
  expect_true(all(feasible$power_actual >= 0.80 - 1e-9))
})

test_that("irrevocable=TRUE is a subset of irrevocable=FALSE", {
  args <- list(p0 = 0.1, p1 = 0.3, alpha = 0.05, power = 0.80, n_max = 20L)
  f_irrev  <- do.call(find_feasible_designs, c(args, irrevocable = TRUE))
  f_all    <- do.call(find_feasible_designs, c(args, irrevocable = FALSE))
  keys <- function(d) paste(d$n1, d$n, d$r1, d$e1, d$r)
  expect_true(all(keys(f_irrev) %in% keys(f_all)))
  expect_lte(nrow(f_irrev), nrow(f_all))
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
  nm  <- admissible$n
  en0 <- admissible$en_null
  en1 <- admissible$en_alt
  for (i in seq_len(m)) {
    dominated <- any(
      seq_len(m) != i &
        nm  <= nm[i]  & en0 <= en0[i] & en1 <= en1[i] &
        (nm  < nm[i]  | en0 < en0[i]  | en1 < en1[i])
    )
    expect_false(dominated, label = sprintf("design %d is dominated", i))
  }
})

test_that("admissible set from synthetic designs has correct Pareto frontier", {
  # All designs have n=20, so 3D Pareto behaves the same as 2D on (en_null, en_alt).
  # D1 dominated by D2 on both EN axes; D3 on frontier (trades EN0 for EN1).
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

test_that("identical (n, en_null, en_alt) triples are all retained", {
  # Neither row strictly dominates the other (no strict improvement on any axis),
  # so both must survive the Pareto filter.
  make_row <- function() {
    data.frame(n1=10, n=20, n2=10, r1=2, e1=8, r=8,
               p0=0.1, p1=0.3,
               alpha_actual=0.04, power_actual=0.82,
               en_null=14, en_alt=12, is_irrevocable=TRUE)
  }
  designs <- rbind(make_row(), make_row())
  adm <- find_admissible_designs(designs)
  expect_equal(nrow(adm), 2L)
})

test_that("identical triples kept; a dominator eliminates all copies", {
  make_row <- function(en0, en1) {
    data.frame(n1=10, n=20, n2=10, r1=2, e1=8, r=8,
               p0=0.1, p1=0.3,
               alpha_actual=0.04, power_actual=0.82,
               en_null=en0, en_alt=en1, is_irrevocable=TRUE)
  }
  designs <- rbind(make_row(15, 12),  # dominated, copy 1
                   make_row(15, 12),  # dominated, copy 2
                   make_row(13, 11))  # dominator
  adm <- find_admissible_designs(designs)
  expect_equal(nrow(adm), 1L)
  expect_equal(adm$en_null, 13)
  expect_equal(adm$en_alt,  11)
})

test_that("3D admissible set retains designs with small n even if worse on EN axes", {
  # D1: small n but worse EN; D2: large n but better EN.
  # In 2D (EN_null, EN_alt) D1 is dominated; in 3D it survives because n[D1] < n[D2].
  make_row <- function(n, en0, en1) {
    data.frame(n1=10, n=n, n2=n-10L, r1=2, e1=8, r=8,
               p0=0.1, p1=0.3,
               alpha_actual=0.04, power_actual=0.82,
               en_null=en0, en_alt=en1, is_irrevocable=TRUE)
  }
  designs <- rbind(make_row(18, 15, 14),  # D1 — small n, worse EN axes
                   make_row(25, 13, 12))  # D2 — large n, better EN axes
  adm <- find_admissible_designs(designs)
  expect_equal(nrow(adm), 2L)  # both survive 3D Pareto
})

# ---------------------------------------------------------------------------
# admissible_designs.R — admissible_designs
# ---------------------------------------------------------------------------

test_that("admissible_designs returns a data.frame with design_type column", {
  result <- admissible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                               power = 0.80, n_max = 25L)
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true("design_type" %in% names(result))
  expect_true(all(result$design_type %in% c("minimax", "optimal",
                                             "minimax, optimal", "")))
})

test_that("exactly one minimax and one optimal label in admissible_designs output", {
  result <- admissible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                               power = 0.80, n_max = 25L)
  mm_rows  <- grepl("minimax", result$design_type)
  opt_rows <- grepl("optimal", result$design_type)
  expect_equal(sum(mm_rows),  1L)
  expect_equal(sum(opt_rows), 1L)
})

test_that("minimax design has the smallest n in the admissible set", {
  result <- admissible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                               power = 0.80, n_max = 25L)
  mm_n <- result$n[grepl("minimax", result$design_type)]
  expect_true(all(mm_n <= result$n))
})

test_that("optimal design has the smallest en_null in the admissible set", {
  result <- admissible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                               power = 0.80, n_max = 25L)
  opt_en <- result$en_null[grepl("optimal", result$design_type)]
  expect_true(all(opt_en <= result$en_null))
})

test_that("admissible_designs with irrevocable=FALSE includes non-irrevocable designs", {
  result <- admissible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                               power = 0.80, n_max = 25L,
                               irrevocable = FALSE)
  expect_gt(nrow(result), 0)
  expect_true(any(!result$is_irrevocable))
})

# ---------------------------------------------------------------------------
# simon = TRUE — Simon two-stage designs (no interim efficacy stop)
# ---------------------------------------------------------------------------

test_that("simon=TRUE returns only designs with e1 = n1 + 1", {
  feasible <- find_feasible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                                    power = 0.80, n_max = 25L, simon = TRUE)
  expect_gt(nrow(feasible), 0)
  expect_true(all(feasible$e1 == feasible$n1 + 1L))
})

test_that("simon=TRUE designs are a subset of all feasible designs", {
  args <- list(p0 = 0.1, p1 = 0.3, alpha = 0.05, power = 0.80, n_max = 20L)
  f_simon <- do.call(find_feasible_designs, c(args, simon = TRUE,
                                              irrevocable = FALSE))
  f_all   <- do.call(find_feasible_designs, c(args, irrevocable = FALSE))
  keys <- function(d) paste(d$n1, d$n, d$r1, d$e1, d$r)
  expect_true(all(keys(f_simon) %in% keys(f_all)))
  expect_lte(nrow(f_simon), nrow(f_all))
})

test_that("admissible_designs with simon=TRUE returns valid output with design_type", {
  result <- admissible_designs(p0 = 0.1, p1 = 0.3, alpha = 0.05,
                               power = 0.80, n_max = 25L, simon = TRUE)
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true("design_type" %in% names(result))
  expect_true(all(result$e1 == result$n1 + 1L))
})
