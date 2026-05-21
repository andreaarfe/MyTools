# twostage

Enumerates all admissible two-stage single-arm clinical trial designs for a
binary endpoint. Supports futility and efficacy stopping at an interim
analysis, with an optional irrevocability constraint ensuring an interim
efficacy declaration cannot be overturned by the final analysis.

## Installation

Install from GitHub with vignettes built locally (requires `knitr` and
`rmarkdown`):

```r
# install.packages("remotes")
remotes::install_github("andreaarfe/MyTools", build_vignettes = TRUE)
```

Without `build_vignettes = TRUE` the package installs but
`vignette("twostage", package = "twostage")` will not find the guide.

## Getting started

```r
library(twostage)
vignette("twostage", package = "twostage")
```

The vignette walks through finding admissible designs, interpreting the
Pareto frontier, and the `irrevocable` and `simon` options.
