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

**Documentation rule:** Whenever a function signature or behaviour changes (parameters, defaults, return value), update both the Roxygen comments in the corresponding `R/*.R` file **and** the generated `man/*.Rd` file. Run `devtools::document()` to regenerate `.Rd` files if R is available; otherwise edit `man/*.Rd` by hand to keep them in sync.

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

**Irrevocability constraint:** `e1 ≥ r` — guarantees that once interim efficacy is declared, no stage-2 outcome can overturn it (worst case: X2 = 0, so X1 + X2 = X1 ≥ e1 ≥ r). Enforced when `irrevocable = TRUE`; default is `FALSE`. Use `TRUE` when the interim efficacy stopping rule is **non-binding** (the trial may continue to stage 2 even after X1 ≥ e1), so that the final decision remains coherent with the interim declaration.

**Simon two-stage designs:** `e1 = n1+1` (no interim efficacy stop). Pass `simon = TRUE` to restrict the search to this class. When `simon = TRUE`, `irrevocable` has no effect — irrevocability is vacuously satisfied because there is never an interim efficacy declaration.

### User-facing documentation

`vignettes/twostage.Rmd` is the primary worked-example guide. It covers: finding admissible designs, interpreting the Pareto frontier, restricting to irrevocable designs (`irrevocable = TRUE`), restricting to Simon two-stage designs (`simon = TRUE`), and evaluating a single design with `evaluate_design()`.

### Module responsibilities

- **`R/design.R`** — thin R wrapper around `evaluate_design_cpp`. Computes `alpha_actual`, `power_actual`, `en_null`, `en_alt`, `is_irrevocable` for a single design. No input validation (caller's responsibility).
- **`R/search.R`** — thin R wrapper around `find_feasible_designs_cpp`. Validates inputs, delegates the five nested loops to C++, emits a message when no designs are found. Returns a `data.frame` with `n` first, `p0`/`p1` dropped, and row names reset. Key options: `irrevocable` (enforce `e1 ≥ r`), `simon` (fix `e1 = n1+1`, no interim efficacy stop).
- **`R/admissibility.R`** — thin R wrapper around `find_admissible_designs_cpp`. Returns the 3D Pareto frontier on `(n, en_null, en_alt)`: a design survives if no other design is weakly better on all three axes with at least one strict improvement. Sorting (by `n` descending) and row-name reset done in R.
- **`R/admissible_designs.R`** — top-level convenience function. Calls `find_feasible_designs` then `find_admissible_designs`, and labels the minimax design (min `n`, tie-break min `en_null`) and optimal design (min `en_null`, tie-break min `n`) via a `design_type` column. Accepts `irrevocable` and `simon` and passes them through.
- **`src/twostage.cpp`** — C++ implementations (`evaluate_design_cpp`, `find_feasible_designs_cpp`, `find_admissible_designs_cpp`). Hot loops live here; binomial PMF cached per `n1`, and survival function `P(X2 ≥ k)` cached per `(n, n1)` to hoist `R::pbinom` out of the innermost loop. `find_admissible_designs_cpp` uses a sort-then-sweep O(m log m) Pareto filter with a 2D staircase rather than the naïve O(m²) double scan.

## Planned work

### Refactor: generalise package scope

Refactor the R package from a two-stage design tool into a general-purpose collection of statistical/clinical-trial utilities. The two-stage design code becomes one module within a broader toolbox; new tools will be added over time. Consider renaming the package and restructuring namespacing accordingly.

### Next feature: OC curves

Add a function to evaluate the operating characteristics of a two-stage design at an arbitrary response probability `p` (not just the design's `p0`/`p1`). This enables operating-characteristic curves across a range of `p` values.

**Outputs of interest:**
- Probability of promotion (final success: `X1 + X2 ≥ r`, or early efficacy: `X1 ≥ e1`)
- Probability of early stopping for futility (`X1 ≤ r1`)
- Probability of early stopping for efficacy (`X1 ≥ e1`)
- Expected final sample size

**Implementation note:** `oc_single()` in `src/twostage.cpp` already computes all four quantities for a given `p` given a precomputed PMF of `X1` and survival table for `X2`. The new R-level function would call `evaluate_design_cpp` (or a thin new C++ helper) in a vectorised loop over a user-supplied vector of `p` values, returning a `data.frame` with one row per `p`. The `build_surv()` helper introduced in the pbinom-caching refactor can be reused directly.

## CRAN compliance

This package targets CRAN. **All changes must be CRAN-compliant.** Key constraints:

- **No custom compiler flags** in `src/Makevars`. Packages must not override the R installation's optimization settings (Writing R Extensions §1.2.1). `-O3`, `-march=native`, and equivalents are prohibited. If you want higher optimization locally, put flags in your personal `~/.R/Makevars` — never in the package's `src/Makevars`.
- **No non-portable C++ extensions.** The bare `restrict` keyword is C99, not C++. Compiler-specific spellings (`__restrict__` on gcc/clang, `__restrict` on MSVC) differ across the toolchains CRAN builds on and must not appear in package source.
- **Must pass `devtools::check()` with 0 errors.** Run `Rscript -e 'devtools::check()'` before committing any change to R or C++ source. Warnings and notes that are pre-existing environment artifacts (locale, qpdf, system R flags) are acceptable; new ones are not.
- **C++ standard.** Stick to C++11 (the Rcpp minimum). Avoid features requiring C++14 or later unless `SystemRequirements` and `CXX_STD` in `DESCRIPTION` are updated and the change is verified to build on all CRAN platforms.

### Known sharp edge

`seq.int(a, b)` in R is **decreasing** when `a > b`, not empty. The continuation region `{r1+1, ..., e1-1}` is empty when `r1 = e1 - 1`; always guard: `if (r1 + 1L <= e1 - 1L) seq.int(...) else integer(0L)`. Omitting this caused ~4× inflation of `alpha_actual` for boundary designs.
