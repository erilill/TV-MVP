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
  T <- nrow(returns)
  p <- ncol(returns)

  residuals <- matrix(NA, nrow = T, ncol = p)

  for (t in 1:T) {
    factors_t <- factors[t, , drop = FALSE]
    loadings_t <- (loadings_list[[t]])

    modeled_returns_t <- factors_t %*% t(loadings_t)
    residuals[t, ] <- returns[t, ] - modeled_returns_t
  }
  return(residuals)
}


#' @import spcov
#' @export
estimate_residual_cov <- function(residuals, tol = 1e-6) {
  p <- ncol(residuals)

  # Compute sample residual covariance matrix regularization
  S_e <- cov(residuals)

  # Initial covariance estimate (diagonal only)
  Sigma_init <- diag(diag(S_e))
  
  lambda <- cv_spcov_lambda(residuals)$optimal_lambda

  # Regularization parameter matrix (lambda is applied to off-diagonal elements)
  Lambda <- matrix(lambda, nrow = p, ncol = p)
  diag(Lambda) <- 0  # Do not penalize diagonal entries

  # Solve using spcov
  result <- spcov(Sigma_init, S_e, lambda = Lambda,
                  step.size=0.001,  trace=0, tol.outer=tol,n.inner.steps = 200,
                  n.outer.steps = 200,thr.inner = 1e-3)

  # Extract estimated covariance matrix
  Sigma_est <- result$Sigma

  return(Sigma_est)
}
#' Compute Time-Varying Covariance Matrices
#'
#' This function computes time-varying covariance matrices based on factor models.
#'
#' @param loadings A matrix of loadings for factors across assets.
#' @param factor_cov A covariance matrix of the factors.
#' @param residual_cov A covariance matrix of the residuals.
#'
#' @return A covariance matrix representing the time-varying covariance of asset returns.
#'
#' @details
#' The function models the covariance matrix of asset returns at each time period \eqn{t}
#' as:
#' \deqn{\Sigma_t = \Lambda_t \Sigma_F \Lambda_t' + \Sigma_{\epsilon}}
#' where \eqn{\Lambda_t} is the loadings matrix, \eqn{\Sigma_F} is the factor covariance matrix,
#' and \eqn{\Sigma_{\epsilon}} is the residual covariance matrix.
#'
#' @examples
#' # Example loadings, factor covariance, and residual covariance
#' loadings <- matrix(runif(25), ncol = 5)
#' factor_cov <- matrix(runif(25), ncol = 5)
#' residual_cov <- diag(0.01, 5)
#'
#' # Compute time-varying covariance
#' time_varying_cov <- compute_time_varying_cov(loadings, factor_cov, residual_cov)
#' print(time_varying_cov)
#'
#' @export
compute_time_varying_cov <- function(loadings, factor_cov, residual_cov) {
  return(loadings %*% factor_cov %*% t(loadings) + residual_cov)
}
