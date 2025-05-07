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
    }

  ),

  private = list(
    # put all the variables here for encapsulation
    # and offer public functions for users to get access them

    # fields that are preserved
    # user cannot set them manually by using set function except data
    important_fields = c("data","iT","ip","important_fields",
                         "optimal_m","IC_values","J_test"),

    data = NULL, # data tibble
    iT = NULL, # integer of the sample size (time)
    ip = NULL, # integer of the number of stocks

    # results
    optimal_m = NULL,
    IC_values = NULL,
    J_test = NULL,

    # user can set this part
    max_m = NULL, # maximum number of factors to consider.
    bandwidth = NULL # bandwidth used in the kernel weighting for the local PCA
  )
)
