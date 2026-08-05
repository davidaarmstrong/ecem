# ecem

Elicited priors for coarsening regimes in coarsened exact matching (CEM).
Companion software for "Elicited Priors for Coarsening Regimes in Coarsened
Exact Matching" (Dery, Alcantara, Armstrong, Stephenson).

CEM requires the analyst to pick a coarsening (a set of cutpoints) for every
covariate before matching. **ecem** lets you instead elicit a *prior* over
those cutpoints -- a range you're genuinely uncertain within, a fixed value
you're not, an exact match, an excluded covariate, or a set of competing
regimes for the same covariate (e.g. a lifecycle vs. a generational theory
of where age should break) -- draw repeatedly from it, match and estimate
on each draw, and pool the results via Rubin's rules. Two diagnostics tell
you whether that pooling is trustworthy, and if it isn't, a fallback
existence test and a set of labeling diagnostics tell you what to report
instead.

## Status

This package was scaffolded by moving working code out of a single script
(`cem_simulation/cem_pipeline.R`) into package form. **None of it has been
run through an R interpreter yet** -- it was written and hand-traced without
R available in that environment. Before relying on it:

```r
devtools::document()   # regenerates NAMESPACE and man/ from the roxygen comments
devtools::load_all()   # loads the package for interactive testing
devtools::check()      # full R CMD check
testthat::test_dir("tests/testthat")
```

Please report anything that errors, or any check/test that fails -- that's
expected at this stage, and the fastest way to find out is to just run it.

## Installation

```r
# not yet on CRAN
devtools::install_local("path/to/ecem")
# or, once pushed to GitHub:
# devtools::install_github("REPLACE-ME/ecem")
```

## A minimal example

```r
library(ecem)

set.seed(1)
pop <- simulate_population(N = 3000, tau0 = 3, heterogeneous = TRUE)

cutpoint_specs <- list(
  age    = list(c(25, 33), c(60, 68)),  # elicited (X_E)
  educ   = c(12, 16),                   # fixed (X_F)
  region = "exact",                     # exact match
  income = NULL,                        # excluded
  noise  = NULL                         # excluded
)

draws  <- run_M_draws(pop, "D", "Y", cutpoint_specs, M = 20)
pooled <- pool_draws(draws)
pooled$tau_bar
```

See `demo(package = "ecem")` for the full workflow, including the
achievable-configuration ceiling / exact enumeration and the competing-
regimes example, and `?ecem` for an overview of every exported function.

## License

MIT
