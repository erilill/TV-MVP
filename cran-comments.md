## Test environments
* local: macOS Sequoia 15.5, R 4.5.0

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Comments
* This is an update to the `TVMVP` package on CRAN. I am the maintainer.
* The package implements a time-varying minimum variance portfolio (TVMVP) model and estimation procedure.
* The package includes documentation, examples, and vignettes.
* Basic tests have been added using `testthat` to validate core functionality (e.g., predict_portfolio()).
* Some examples are wrapped in \donttest{} because they involve time-consuming portfolio optimization routines using. These are provided for reproducibility but would exceed the recommended runtime for CRAN checks.