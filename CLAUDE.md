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
# or, if devtools is unavailable (e.g. remote web environment):
Rscript -e 'roxygen2::roxygenise(".")'
```

**Documentation rule:** Whenever a function signature or behaviour changes (parameters, defaults, return value), update both the Roxygen comments in the corresponding `R/*.R` file **and** the generated `man/*.Rd` file. Regenerate `.Rd` files with `devtools::document()` (preferred) or `roxygen2::roxygenise(".")` (available via apt as `r-cran-roxygen2`); never edit `man/*.Rd` by hand.

### Remote web environment (Claude Code on the web)

In this environment R 4.3.3 is available via `apt`, but CRAN is unreachable so `devtools` cannot be installed. Use these equivalents:

| devtools command | alternative |
|---|---|
| `devtools::test()` | `pkgload::load_all("."); testthat::test_local(".")` |
| `devtools::load_all()` | `pkgload::load_all(".")` |
| `devtools::document()` | `roxygen2::roxygenise(".")` |
| `devtools::check()` | `R CMD build . && R CMD check --as-cran MyTools_*.tar.gz` |

`r-cran-rcpp`, `r-cran-testthat`, `r-cran-pkgload`, and `r-cran-roxygen2` are all apt-packaged and installable with `sudo apt-get install -y`.

**After editing `src/twostage.cpp`** (regenerates R/RcppExports.R and src/RcppExports.cpp):
```bash
Rscript -e 'Rcpp::compileAttributes()'
```

## Architecture

`MyTools` is a collection of statistical and clinical-trial utilities, organised as a set of named modules. Each module is identified by a function-name prefix (e.g. `twostage_*`) and a set of files prefixed by the same name in `R/`, `man/`, and the vignettes index. Future modules will follow the same pattern; shared C++ utilities will be extracted into their own source files when a second module needs them.

### Module: twostage

The `twostage` module enumerates all admissible two-stage single-arm clinical trial designs for a binary endpoint. The pipeline is: **search → filter → Pareto-optimise**.

### Core design parameters

Each design is defined by `(n1, n, r1, e1, r)`:
- `n1` / `n`: stage-1 and total sample sizes
- `r1`: futility boundary — stop for futility if `X1 ≤ r1` (`-1` = no futility stop)
- `e1`: interim efficacy boundary — declare success if `X1 ≥ e1` (`n1+1` = no interim stop)
- `r`: final success threshold — success if `X1 + X2 ≥ r`

**Non-binding constraint:** `e1 ≥ r` — guarantees that once interim efficacy is declared, no stage-2 outcome can overturn it (worst case: X2 = 0, so X1 + X2 = X1 ≥ e1 ≥ r). Enforced when `non_binding = TRUE`; default is `FALSE`. Use `TRUE` when the interim efficacy stopping rule is **non-binding** (the trial may continue to stage 2 even after X1 ≥ e1), so that the final decision remains coherent with the interim declaration.

**Simon two-stage designs:** `e1 = n1+1` (no interim efficacy stop). Pass `simon = TRUE` to restrict the search to this class. When `simon = TRUE`, `non_binding` has no effect — the `e1 ≥ r` constraint is vacuously satisfied because there is never an interim efficacy declaration.

### User-facing documentation

`vignettes/twostage.Rmd` is the primary worked-example guide for the `twostage` module. It covers: finding admissible designs, interpreting the Pareto frontier, restricting to non-binding designs (`non_binding = TRUE`), restricting to Simon two-stage designs (`simon = TRUE`), and evaluating a design at one or more response probabilities with `twostage_evaluate_design()`.

### Module responsibilities

- **`R/twostage-design.R`** — thin R wrapper around `evaluate_design_cpp` exporting `twostage_evaluate_design()`. Evaluates a fixed design across a vector of response probabilities `p`, returning a `data.frame` with columns `p`, `p_success`, `p_futility`, `p_efficacy1`, and `en`. Scalar `p` yields a one-row frame.
- **`R/twostage-search.R`** — internal (not exported) R wrapper around `find_feasible_designs_cpp`. Validates inputs, delegates the five nested loops to C++, emits a message when no designs are found. Returns a `data.frame` with `n` first, `p0`/`p1` dropped, and row names reset. Key options: `non_binding` (enforce `e1 ≥ r`), `simon` (fix `e1 = n1+1`, no interim efficacy stop). Called only by `twostage_admissible_designs()`; tests reach it via `MyTools:::twostage_find_feasible_designs()`.
- **`R/twostage-admissibility.R`** — internal (not exported) R wrapper around `find_admissible_designs_cpp`. Returns the 3D Pareto frontier on `(n, en_null, en_alt)`: a design survives if no other design is weakly better on all three axes with at least one strict improvement. Sorting (by `n` descending) and row-name reset done in R. Called only by `twostage_admissible_designs()`; tests reach it via `MyTools:::twostage_find_admissible_designs()`.
- **`R/twostage-admissible_designs.R`** — top-level convenience function `twostage_admissible_designs()`. Calls `twostage_find_feasible_designs` then `twostage_find_admissible_designs`, and labels four designs of interest (minimax, optimal null, optimal alt, optimal mean) via a `design_type` column. Accepts `non_binding` and `simon` and passes them through.
- **`src/twostage.cpp`** — C++ implementations (`evaluate_design_cpp`, `find_feasible_designs_cpp`, `find_admissible_designs_cpp`). Hot loops live here; binomial PMF cached per `n1`, and survival function `P(X2 ≥ k)` cached per `(n, n1)` to hoist `R::pbinom` out of the innermost loop. `find_admissible_designs_cpp` uses a sort-then-sweep O(m log m) Pareto filter with a 2D staircase rather than the naïve O(m²) double scan.

## CRAN compliance

This package targets CRAN. **All changes must be CRAN-compliant.** Key constraints:

- **No custom compiler flags** in `src/Makevars`. Packages must not override the R installation's optimization settings (Writing R Extensions §1.2.1). `-O3`, `-march=native`, and equivalents are prohibited. If you want higher optimization locally, put flags in your personal `~/.R/Makevars` — never in the package's `src/Makevars`.
- **No non-portable C++ extensions.** The bare `restrict` keyword is C99, not C++. Compiler-specific spellings (`__restrict__` on gcc/clang, `__restrict` on MSVC) differ across the toolchains CRAN builds on and must not appear in package source.
- **Must pass `devtools::check()` with 0 errors.** Run `Rscript -e 'devtools::check()'` before committing any change to R or C++ source. Warnings and notes that are pre-existing environment artifacts (locale, qpdf, system R flags) are acceptable; new ones are not.
- **C++ standard.** Stick to C++11 (the Rcpp minimum). Avoid features requiring C++14 or later unless `SystemRequirements` and `CXX_STD` in `DESCRIPTION` are updated and the change is verified to build on all CRAN platforms.

### Known sharp edge

`seq.int(a, b)` in R is **decreasing** when `a > b`, not empty. The continuation region `{r1+1, ..., e1-1}` is empty when `r1 = e1 - 1`; always guard: `if (r1 + 1L <= e1 - 1L) seq.int(...) else integer(0L)`. Omitting this caused ~4× inflation of `alpha_actual` for boundary designs.
