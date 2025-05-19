## Test environments
* local: macOS Sequoia 15.5, R 4.5.0

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:

* checking for future file timestamps ... NOTE
  unable to verify current time

This NOTE is caused by limitations in the system clock resolution and is not related to the package itself. It does not affect reproducibility or functionality.

## Additional information
* This is the initial submission of the `TVMVP` package to CRAN.
* The package implements a time-varying minimum variance portfolio (TVMVP) model and estimation procedure as introduced in [Your Reference or Paper if applicable].
* The package includes documentation and examples.
* Basic tests have been added using `testthat` to validate core functionality (e.g., predict_portfolio()).
