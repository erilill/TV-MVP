## In this file, we define the R6 type class of the package

#' Time Varying Minimum Variance Portfolio (TVMVP) Class
#'
#' This class implements a time-varying mean-variance portfolio model.
#'
#' @export
TVMVP <- R6::R6Class(
  "TVMVP",
  
  public = list(
    
    initialize = function(data = NULL) {
      if(is_tibble(data)){
        self$data <- data
      } esle{
      }
    }
    
  ),
  
  private = list(
    # put all the variables here for encapsulation
    # and offer public functions for users to get access them
    
    data <- NULL
    
  )
)