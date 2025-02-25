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
      if(tibble::is_tibble(data)){
        self$data <- data
        cli::cli_alert_success("Tibble data set successfully loaded.")
        # can also give the information of the data set
        # for example, the size in MB, row and column numbers etc.
      } esle{
        cli_alert_info("The data set is not tibble!")
      }
    }
    
  ),
  
  private = list(
    # put all the variables here for encapsulation
    # and offer public functions for users to get access them
    
    data <- NULL
    
  )
)