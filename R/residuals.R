#' Estimate Residuals from Factor Model
#'
#' This function estimates the residuals of asset returns after removing the effect
#' of factor-driven returns.
#'
#' @param factors_list A list where each element is a matrix representing factor data
#' for a specific time period.
#' @param loadings_list A list where each element is a matrix of loadings corresponding
#' to the factors for each time period.
#' @param returns A matrix of asset returns with rows representing time periods and
#' columns representing assets.
#'
#' @return A matrix of residuals where each row corresponds to a time period and each
#' column corresponds to an asset.
#'
#' @details
#' For each time period \eqn{t}, the function models the asset returns as:
#' \deqn{R_t = F_t \Lambda_t + \epsilon_t}
#' where \eqn{R_t} is the vector of asset returns, \eqn{F_t} is the factor matrix,
#' \eqn{\Lambda_t} is the loadings matrix, and \eqn{\epsilon_t} represents the residuals.
#'
#' The residuals are computed as the difference between actual returns and the modeled
#' returns.
#'
#' @examples
#' # Example factors and loadings
#' factors_list <- list(
#'   matrix(rnorm(100), ncol = 5),
#'   matrix(rnorm(100), ncol = 5)
#' )
#' loadings_list <- list(
#'   matrix(runif(25), ncol = 5),
#'   matrix(runif(25), ncol = 5)
#' )
#'
#' # Simulate returns matrix
#' returns <- matrix(rnorm(200 * 5, mean = 0.01, sd = 0.02), ncol = 5)
#'
#' # Estimate residuals
#' resid <- residuals(factors_list, loadings_list, returns)
#' head(resid)
#'
#' @export
residuals <- function(factors, loadings_list, returns) {
  iT <- nrow(returns)
  ip <- ncol(returns)
  
  residuals <- matrix(NA, nrow = iT, ncol = ip)
  
  for (t in 1:iT) {
    factors_t <- factors[t, , drop = FALSE]
    loadings_t <- (loadings_list[[t]])
    
    modeled_returns_t <- factors_t %*% t(loadings_t)
    residuals[t, ] <- returns[t, ] - modeled_returns_t
  }
  return(residuals)
}
