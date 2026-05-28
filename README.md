# MyTools

A collection of R utilities for statistical and clinical-trial methodology.
The current release provides the **twostage** module (functions prefixed
`twostage_`), which enumerates admissible two-stage single-arm trial designs
for a binary endpoint, with optional non-binding and Simon-design
constraints.

## Installation

Install from GitHub with vignettes built locally (requires `knitr` and
`rmarkdown`):

```r
# install.packages("remotes")
remotes::install_github("andreaarfe/MyTools", build_vignettes = TRUE)
```

Without `build_vignettes = TRUE` the package installs but
`vignette("twostage", package = "MyTools")` will not find the guide.

## Getting started

```r
library(MyTools)
vignette("twostage", package = "MyTools")
```

The vignette walks through finding admissible designs, interpreting the
Pareto frontier, and the `non_binding` and `simon` options.

## Modules

- **twostage** (`twostage_*` functions): admissible two-stage single-arm
  trial designs for a binary endpoint.
