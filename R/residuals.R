#' Estimate Residuals from Factor Model
#'
#' This function estimates the residuals of asset returns after removing the effect
#' of factor-driven returns.
#'
#' @param factors A matrix containing the step-ahead-factors of from the \code{localPCA} function.
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
#' where \eqn{R_t} is the vector of asset returns, \eqn{F_t} is the t'th row of the factor matrix,
#' \eqn{\Lambda_t} is the loadings matrix, and \eqn{\epsilon_t} represents the residuals.
#'
#' The residuals are computed as the difference between actual returns and the modeled
#' returns.
#'
#' 
#' @keywords internal
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
