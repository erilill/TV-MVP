## In this file, we define the R6 type class of the package

#' Time Varying Minimum Variance Portfolio (TVMVP) Class
#'
#' @description
#' This class implements a time-varying minimum variance portfolio using locally
#' smoothed principal component analysis (PCA) to estimate the time-dependent
#' covariance matrix.
#' 
#' This class provides a flexible interface to:
#' \itemize{
#'    \item Set return data (\code{$set_data()})
#'    \item Determine optimal number of factors ({\code{$\link{determine_factors}()}})
#'    \item Conduct test of constant factor loadings (\code{$\link{hyptest}()})
#'    \item Time-dependent covariance estimation (\code{$\link{time_varying_cov}()})
#'    \item Portfolio optimization (\code{$\link{predict_portfolio}()})
#'    \item Rolling window evaluation (\code{$\link{rolling_time_varying_mvp}()})
#'    \item Extract cached results (\code{$get_optimal_m()}, \code{$get_IC_values()}, 
#'    \code{$get_bootstrap()})
#' }
#' 
#' Looking for package description? See \link{TVMVP-package}.
#'
#' @section Usage:
#' \preformatted{
#' # Initial object of class TVMVP
#' tv <- TVMVP$new()
#' 
#' # Set data
#' tv$set_data(returns) # Returns must be T times p matrix
#' 
#' # Determine number of factors
#' tv$determine_factors(max_m=10)
#' 
#' # Test for constant loadings
#' tv$hyptest()
#' 
#' # Estimate time-dependent covariance matrix
#' cov <- tv$time_varying_cov()
#' 
#' # Evaluate TVMVP performance on historical data
#' mvp_results <- tv$rolling_time_varying_mvp(
#'                  initial_window = 60,
#'                  rebal_period   = 5,
#'                  max_factors    = 10,
#'                  return_type    = "daily")
#'                  
#' # Make out-of-sample prediction and compute weights
#' predictions <- tv$predict_portfolio(horizon=5,
#'                                    min_return = 0.01,
#'                                    max_SR = TRUE)
#' # Extract weights
#' predictions$getWeights("MVP")
#' }
#' #'
#' 
#' @section Arguments: 
#' \describe{
#'    \item{data}{\eqn{T × p} (time periods by assets) matrix of returns.}
#'    \item{bandwidth}{Numerical. Bandwidth parameter used for local smoothing in the local PCA}
#'    \item{max_m}{Integer. Maximum number of factors to be tested when determining the optimal number of factors.}
#'    \item{optimal_m}{The optimal number of factors to use in covariance estimation.}
#'  }
#' 
#'
#' @section Methods:
#' \describe{
#'    \item{\code{$new(data = NULL)}}{Initialize object of class TVMVP. Optionally pass returns matrix.}
#'    \item{\code{$set_data(data)}}{Set the data. Must be \eqn{T × p} (time periods by assets) matrix.}
#'    \item{\code{$get_data()}}{Get the data.}
#'    \item{\code{$set()}}{Manually set arguments of the object.}
#'    \item{\code{$\link{determine_factors}()}}{Determines optimal number of factors based on BIC-type information criterion.}
#'    \item{\code{$get_optimal_m{}}}{Prints optimal number of factors, \code{optimal_m.}}
#'    \item{\code{$get_IC_values()}}{Prints IC-values for the different number of factors tested using \code{\link{determine_factors}}.}
#'    \item{\code{$\link[=hyptest1]{hyptest()}}}{Hypothesis test of constant loadings.}
#'    \item{\code{$get_bootstrap()}}{Prints bootstrap test statistics from the hypothesis test.}
#'    \item{\code{$\link{predict_portfolio}()}}{Optimizes portfolio weights for out-of-sample prediction of portfolio performance.}
#'    \item{\code{$\link{rolling_time_varying_mvp}()}}{Evaluates MVP performance in a rolling window framework.}
#'    \item{\code{$\link{time_varying_cov}()}}{Estimates the time-varying covariance matrix.}
#'    \item{\code{$\link{silverman}()}}{Silverman's rule of thumb bandwidth formula.}
#' }
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
