## Test environments
* local: macOS Sequoia 15.5, R 4.5.0

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Comments
* This is the initial submission of the `TVMVP` package to CRAN. I am the maintainer.
* The package implements a time-varying minimum variance portfolio (TVMVP) model and estimation procedure as introduced in [Your Reference or Paper if applicable].
* The package includes documentation, examples, and vignettes.
* Basic tests have been added using `testthat` to validate core functionality (e.g., predict_portfolio()).
* The NOTE about "Author field differs from that derived from Authors@R" is expected. All author metadata is provided via Authors@R with ORCID and affiliation in comment = ..., following current R packaging guidelines.
* Some examples are wrapped in \dontrun{} because they involve time-consuming portfolio optimization routines using. These are provided for reproducibility but would exceed the recommended runtime for CRAN checks.