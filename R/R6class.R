## In this file, we define the R6 type class of the package

#' Time Varying Minimum Variance Portfolio (TVMVP) Class
#'
#' @description
#' This class implements a time-varying mean-variance portfolio model.
#'
#' @section Usage:
#' ```
#' ```
#' @section Arguments:
#' * `data`: A numeric matrix for the data matrix.
#'
#'
#' @section Methods:
#' - `$new(data = NULL)`: Initialize the class.
#'
#' - `$set_data(data)`: Set the data.
#'
#' - `$get_data()`: Get the data.
#'
#' @name TVMVP
#'
#' @docType class
NULL

#' @import R6
#' @import cli
#' @import prettyunits
#' @export
TVMVP <- R6::R6Class(
  "TVMVP",
  cloneable = FALSE,
  public = list(


    initialize = function(data = NULL) {
      if(!is.null(data)) self$set(data=data)
    },

    determine_factors = function(max_R, bandwidth) {
      p_determine_factors(self, private, max_R, bandwidth)
    }
  ),

  private = list(
    # put all the variables here for encapsulation
    # and offer public functions for users to get access them

    data = NULL, # data tibble
    iT = NULL, # integer of the sample size (time)
    ip = NULL, # integer of the number of stocks

    # user can set this part
    max_m = NULL, # maximum number of factors to consider.
    bandwidth = NULL # bandwidth used in the kernel weighting for the local PCA
  )
)
