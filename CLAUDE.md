# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Run all tests** (from repo root):
```bash
Rscript -e 'devtools::test()'
```

**Run a single test by name:**
```bash
Rscript -e 'devtools::test(filter = "empty continuation region")'
```

**Run R CMD check:**
```bash
Rscript -e 'devtools::check()'
```

**Load all modules interactively** (from repo root):
```r
devtools::load_all()
```

**Run the end-to-end example** (after installing or loading the package):
```bash
Rscript examples/example_design.R
```

**Build vignettes:**
```bash
Rscript -e 'devtools::build_vignettes()'
```

**After editing Roxygen comments in any `R/*.R` file** (regenerates `man/*.Rd`):
```bash
Rscript -e 'devtools::document()'
```

**After editing `src/twostage.cpp`** (regenerates R/RcppExports.R and src/RcppExports.cpp):
```bash
Rscript -e 'Rcpp::compileAttributes()'
```

## Architecture

The package enumerates all admissible two-stage single-arm clinical trial designs for a binary endpoint. The pipeline is: **search → filter → Pareto-optimise**.

### Core design parameters

Each design is defined by `(n1, n, r1, e1, r)`:
- `n1` / `n`: stage-1 and total sample sizes
- `r1`: futility boundary — stop for futility if `X1 ≤ r1` (`-1` = no futility stop)
- `e1`: interim efficacy boundary — declare success if `X1 ≥ e1` (`n1+1` = no interim stop)
- `r`: final success threshold — success if `X1 + X2 ≥ r`

**Irrevocability constraint:** `e1 ≥ r` — guarantees that once interim efficacy is declared, no stage-2 outcome can overturn it (worst case: X2 = 0, so X1 + X2 = X1 ≥ e1 ≥ r). Enforced when `irrevocable = TRUE` (default).

**Simon two-stage designs:** `e1 = n1+1` (no interim efficacy stop). Pass `simon = TRUE` to restrict the search to this class. When `simon = TRUE`, `irrevocable` has no effect — irrevocability is vacuously satisfied because there is never an interim efficacy declaration.

### User-facing documentation

`vignettes/twostage.Rmd` is the primary worked-example guide. It covers: finding admissible designs, interpreting the Pareto frontier, searching without the irrevocability constraint (`irrevocable = FALSE`), restricting to Simon two-stage designs (`simon = TRUE`), and evaluating a single design with `evaluate_design()`.

### Module responsibilities

- **`R/distributions.R`** — stateless wrappers: `binom_pmf`, `binom_cdf`, `binom_upper`. Use `pbinom(..., lower.tail = FALSE)` for numerical stability.
- **`R/design.R`** — thin R wrapper around `evaluate_design_cpp`. Computes `alpha_actual`, `power_actual`, `en_null`, `en_alt`, `is_irrevocable` for a single design. No input validation (caller's responsibility).
- **`R/search.R`** — thin R wrapper around `find_feasible_designs_cpp`. Validates inputs, delegates the five nested loops to C++, emits a message when no designs are found. Returns a `data.frame`. Key options: `irrevocable` (enforce `e1 ≥ r`), `simon` (fix `e1 = n1+1`, no interim efficacy stop).
- **`R/admissibility.R`** — thin R wrapper around `find_admissible_designs_cpp`. Returns the 3D Pareto frontier on `(n, en_null, en_alt)`: a design survives if no other design is weakly better on all three axes with at least one strict improvement. Sorting (by `en_null`) done in R.
- **`R/admissible_designs.R`** — top-level convenience function. Calls `find_feasible_designs` then `find_admissible_designs`, and labels the minimax design (min `n`, tie-break min `en_null`) and optimal design (min `en_null`, tie-break min `n`) via a `design_type` column. Accepts `irrevocable` and `simon` and passes them through.
- **`src/twostage.cpp`** — C++ implementations (`evaluate_design_cpp`, `find_feasible_designs_cpp`, `find_admissible_designs_cpp`). Hot loops live here; binomial PMF cached per `n1` inside the search loop.

### Known sharp edge

`seq.int(a, b)` in R is **decreasing** when `a > b`, not empty. The continuation region `{r1+1, ..., e1-1}` is empty when `r1 = e1 - 1`; always guard: `if (r1 + 1L <= e1 - 1L) seq.int(...) else integer(0L)`. Omitting this caused ~4× inflation of `alpha_actual` for boundary designs.
