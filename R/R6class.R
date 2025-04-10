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
                                    cli::cli_text("  - GMV: Use object$GMV")
                                    if (!is.null(self$max_SR)) {
                                      cli::cli_text("  - Maximum Sharpe Ratio Portfolio: Use object$max_SR")
                                    }
                                    if (!is.null(self$MinVarWithReturnConstraint)) {
                                      cli::cli_text("  - Minimum Variance Portfolio with Return Constraint: Use object$MinVarWithReturnConstraint")
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
#' @import R6
RollingWindow <- R6::R6Class(
  "RollingWindow",
  public = list(
    summary = NULL,
    TVMVP   = NULL,
    Equal  = NULL,
    
    initialize = function(summary, TVMVP, Equal) {
      self$summary <- summary
      self$TVMVP   <- TVMVP
      self$Equal   <- Equal
    },
    
    print = function(...) {
      # Header
      cli::cli_h1("Rolling Window Portfolio Analysis")
      cli::cli_rule()
      
      # Print summary
      cli::cli_h2("Summary Metrics")
      df <- self$summary
      
      # Here, if you want, you can rename the 'Method' column so it prints
      # nicely, or leave it as is. For example:
      # df$Method <- with(df, ifelse(Method == "Time-Varying MVP", 
      #                              "Time-Varying MVP Portfolio",
      #                       ifelse(Method == "Equal Weight", 
      #                              "Equal-Weighted Portfolio", 
      #                              Method)))
      
      print(df, row.names = FALSE)
      
      cli::cli_rule()
      cli::cli_h2("Detailed Components")
      cli::cli_text("The detailed portfolio outputs are stored in the following elements:")
      cli::cli_text("  - Time-Varying MVP: Access via `$TVMVP`")
      cli::cli_text("  - Equal Weight: Access via `$Equal`")
      
      invisible(self)
    },
    
    getWeights = function(method = c("TVMVP", "Equal")) {
      method <- match.arg(method)
      if (method == "TVMVP") {
        return(self$TVMVP$weights)
      } else {
        return(self$Equal$weights)
      }
    }
  )
)
