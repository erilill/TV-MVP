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
residuals <- function(factors, loadings_list, returns){
  T <- nrow(returns)
  p <- ncol(returns)
  m <- dim(factors)[2]

  residuals <- matrix(NA, nrow = T, ncol = p)
  for (t in 1:T) {
    modeled_returns_t <- factors[t,] %*% t(loadings_list[[t]])
    residuals[t, ] <- returns[t, ] - modeled_returns_t
  }
  return(residuals)
}
#' Estimate Residual Covariance Matrix with Lasso Penalization
#'
#' This function estimates the residual covariance matrix from asset returns using
#' Lasso penalization to induce sparsity.
#'
#' @param residuals A matrix of residuals where rows represent time periods and columns
#' represent assets.
#' @param lambda A numeric value specifying the penalty parameter for Lasso regression.
#' Defaults to \code{0.1}.
#'
#' @return A sparse covariance matrix estimated from the residuals.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the sample covariance matrix of residuals.
#'   \item Initializes a sparse covariance matrix with zeros.
#'   \item Applies Lasso regression row-wise to estimate the sparse structure of the
#'   covariance matrix.
#'   \item Symmetrizes the covariance matrix and retains the original diagonal elements.
#' }
#'
#' @examples
#' # Simulate residuals matrix
#' residuals_matrix <- matrix(rnorm(200 * 5, mean = 0, sd = 0.01), ncol = 5)
#'
#' # Estimate residual covariance with Lasso
#' residual_cov <- estimate_residual_cov(residuals_matrix, lambda = 0.1)
#' print(residual_cov)
#'
#' @import glmnet
#' @export
estimate_residual_cov <- function(residuals, lambda = 0.1) {
  p <- ncol(residuals)

  # Compute the sample covariance matrix
  S <- cov(residuals)

  # Initialize the sparse covariance matrix
  sparse_cov <- matrix(0, nrow = p, ncol = p)

  # Apply Lasso regression row-wise
  for (i in 1:p) {
    response <- S[i, -i]
    predictors <- S[-i, -i]

    # Fit Lasso regression without intercept
    fit <- glmnet::glmnet(predictors, response, alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE)

    # Extract coefficients and assign to sparse covariance matrix
    coef_i <- as.vector(coef(fit, s = lambda))[-1]  # Exclude intercept
    sparse_cov[i, -i] <- coef_i
  }

  # Symmetrize the covariance matrix
  sparse_cov <- (sparse_cov + t(sparse_cov)) / 2

  # Assign original diagonal elements
  diag(sparse_cov) <- diag(S)

  return(sparse_cov)
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
