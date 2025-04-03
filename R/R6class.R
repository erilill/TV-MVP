## In this file, we define the R6 type class of the package

#' Time Varying Minimum Variance Portfolio (TVMVP) Class
#'
#' This class implements a time-varying mean-variance portfolio model.
#'
#' @import R6
#' @export
TVMVP <- R6::R6Class(
  "TVMVP",

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

    # set functions

    set_data = function(data){
      # for set_data function, the argument data should not be missing
      # this is why I use missing(data) for checking without setting any default value of it
      # note that this is different from initialize function where data can be missing with default value NULL
      if(missing(data)){
        cli::cli_alert_warning("You forgot input the data!")
      } else{
        private$data <- data
        private$iT <- nrow(data)
        private$ip <- ncol(data)

        #tmp_size <- get_object_size(private$data)
        tmp_size <- prettyunits::pretty_bytes(object.size(private$data))
        cli::cli_alert_info("Tibble data set {tmp_size} with {private$iT} rows and {private$ip} columns successfully assigned.")
      }

      invisible(self)  # Enables chaining
    },

    # get functions

    get_data = function(){
      return(private$data)
      # make a copy of the data and return it
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
#' @import R6
PortfolioPredictions <- R6Class("PortfolioPredictions",
                                public = list(
                                  summary = NULL,
                                  GMV = NULL,
                                  max_SR = NULL,
                                  MinVarWithReturnConstraint = NULL,
                                  
                                  initialize = function(summary, GMV, max_SR = NULL, MinVarWithReturnConstraint = NULL) {
                                    self$summary <- summary
                                    self$GMV <- GMV
                                    self$max_SR <- max_SR
                                    self$MinVarWithReturnConstraint <- MinVarWithReturnConstraint
                                  },
                                  
                                  print = function(...) {
                                    cli::cli_h1("Portfolio Optimization Predictions")
                                    cli::cli_rule()
                                    cli::cli_h2("Summary Metrics")
                                    df <- self$summary
                                    df$Method <- with(df, ifelse(Method == "GMV", 
                                                                 "Minimum Variance Portfolio", 
                                                                 ifelse(Method == "max_SR", 
                                                                        "Maximum SR Portfolio", 
                                                                        ifelse(Method == "MinVarWithReturnConstraint", 
                                                                               "Return-Constrained Portfolio", Method))))
                                    print(df, row.names = FALSE)
                                    cli::cli_rule()
                                    cli::cli_h2("Detailed Components")
                                    cli::cli_text("The detailed portfolio outputs are stored in the following elements:")
                                    cli::cli_text("  • GMV: Use object$GMV")
                                    if (!is.null(self$max_SR)) {
                                      cli::cli_text("  • Maximum Sharpe Ratio Portfolio: Use object$max_SR")
                                    }
                                    if (!is.null(self$MinVarWithReturnConstraint)) {
                                      cli::cli_text("  • Minimum Variance Portfolio with Return Constraint: Use object$MinVarWithReturnConstraint")
                                    }
                                    invisible(self)
                                  },
                                  
                                  getWeights = function(method = "GMV") {
                                    switch(method,
                                           GMV = self$GMV$weights,
                                           max_SR = {
                                             if (!is.null(self$max_SR)) {
                                               self$max_SR$weights
                                             } else {
                                               cli::cli_alert_danger("max_SR portfolio not available!")
                                               NULL
                                             }
                                           },
                                           MinVarWithReturnConstraint = {
                                             if (!is.null(self$MinVarWithReturnConstraint)) {
                                               self$MinVarWithReturnConstraint$weights
                                             } else {
                                               cli::cli_alert_danger("MinVarWithReturnConstraint portfolio not available!")
                                               NULL
                                             }
                                           },
                                           stop("Method not found")
                                    )
                                  }
                                )
)
