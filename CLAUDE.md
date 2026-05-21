# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Run all tests** (from repo root):
```bash
Rscript -e 'testthat::test_file("tests/test_design.R")'
```

**Run a single test by name:**
```bash
Rscript -e 'testthat::test_file("tests/test_design.R", filter = "empty continuation region")'
```

**Run the end-to-end example:**
```bash
Rscript examples/example_design.R
```

**Source all modules interactively** (from repo root):
```r
source("R/distributions.R")
source("R/design.R")
source("R/search.R")
source("R/admissibility.R")
```

## Architecture

The package enumerates all admissible two-stage single-arm clinical trial designs for a binary endpoint. The pipeline is: **search → filter → Pareto-optimise**.

### Core design parameters

Each design is defined by `(n1, n, r1, e1, r)`:
- `n1` / `n`: stage-1 and total sample sizes
- `r1`: futility boundary — stop for futility if `X1 ≤ r1` (`-1` = no futility stop)
- `e1`: interim efficacy boundary — declare success if `X1 ≥ e1` (`n1+1` = no interim stop)
- `r`: final success threshold — success if `X1 + X2 ≥ r`

**Irrevocability constraint:** `e1 ≥ r` — guarantees that once interim efficacy is declared, no stage-2 outcome can overturn it (worst case: X2 = 0, so X1 + X2 = X1 ≥ e1 ≥ r).

### Module responsibilities

- **`R/distributions.R`** — stateless wrappers: `binom_pmf`, `binom_cdf`, `binom_upper`. Use `pbinom(..., lower.tail = FALSE)` for numerical stability.
- **`R/design.R`** — `evaluate_design(n1, n, r1, e1, r, p0, p1)` computes `alpha_actual`, `power_actual`, `en_null`, `en_alt`, `is_irrevocable`. Returns a named list; no validation (caller filters).
- **`R/search.R`** — `find_feasible_designs(p0, p1, alpha, power, n_max)` loops over all `(n, n1, r, e1, r1)` combinations. Key pruning: skip `r > n1` because `e1 ≤ n1 < r` would violate irrevocability for all `e1`. Returns a `data.frame`.
- **`R/admissibility.R`** — `find_admissible_designs(designs)` returns the Pareto frontier on `(en_null, en_alt)`: a design survives if no other design is ≤ on both axes with at least one strict inequality.

### Known sharp edge

`seq.int(a, b)` in R is **decreasing** when `a > b`, not empty. The continuation region `{r1+1, ..., e1-1}` is empty when `r1 = e1 - 1`; always guard: `if (r1 + 1L <= e1 - 1L) seq.int(...) else integer(0L)`. Omitting this caused ~4× inflation of `alpha_actual` for boundary designs.
