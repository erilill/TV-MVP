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
#' * `data`: A tibble matrix for the data matrix.
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
#' @import tibble
#' @export
TVMVP <- R6::R6Class(
  "TVMVP",
  cloneable = FALSE,
  public = list(

    
    initialize = function(data = NULL) {
      if(!is.null(data)){
        if(tibble::is_tibble(data)){
          private$data <- data
          private$iT <- nrow(data)
          private$ip <- ncol(data)

          #tmp_size <- get_object_size(private$data)
          tmp_size <- prettyunits::pretty_bytes(object.size(private$data))
          cli::cli_alert_info("Tibble data set {tmp_size} with {private$iT} rows and {private$ip} columns successfully assigned.")
        } else{
          cli_alert_info("The data set is not tibble! The data is empty now.")
        }
      }
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
    ip = NULL # integer of the number of stocks

  )
)