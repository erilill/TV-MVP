#' Forecast Local Factor Model
#'
#' This function implements a local factor model to forecast asset returns
#' using a time-varying covariance structure estimated through kernel-weighted PCA.
#'
#' @param returns A matrix of asset returns, where rows represent time periods (T) and columns represent assets (p).
#' @param W An integer specifying the rolling window size used for estimation.
#' @param m An integer specifying the maximum number of factors to use in the PCA.
#' @param bandwidth A numeric value for the kernel bandwidth, used to calculate the boundary kernel weights.
#' @param kernel_func A kernel function used for weighting. The default is \code{\link{epanechnikov_kernel}}.
#' @param lambda A numeric penalization parameter used in the estimation of residual covariance matrices. Default is 0.1.
#'
#' @details The function forecasts returns using the following steps:
#' \enumerate{
#'   \item Computes kernel weights for past observations.
#'   \item Applies weighted PCA to extract factors and loadings.
#'   \item Estimates factor and residual covariance matrices.
#'   \item Combines factor and residual covariances to compute the time-varying covariance matrix.
#'   \item Computes optimal portfolio weights based on the estimated covariance matrix.
#' }
#'
#' The forecasts are generated for time periods from \code{(W + 1)} to \code{T},
#' as earlier periods are used for model initialization.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{forecasts}}{A T x p matrix of forecasted returns, where \code{NA} values are present for periods \code{1:W}.}
#'   \item{\code{est_covariances}}{A list of length T, where each element is a p x p estimated covariance matrix.}
#'   \item{\code{residual_covariances}}{A list of length T, where each element is a p x p estimated residual covariance matrix.}
#'   \item{\code{optimal_weights}}{A list of length T, where each element is a p-dimensional vector of portfolio weights.}
#' }
#'
#' @seealso \code{\link{epanechnikov_kernel}}, \code{\link{estimate_residual_cov}}, \code{\link{compute_optimal_weights}}
#'
#' @examples
#' set.seed(123)
#' T <- 200  # Number of time periods
#' p <- 50   # Number of assets
#' returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
#' bandwidth <- 0.1
#' W <- 50
#' m <- 3
#'
#' results <- forecast_local_factor_model(
#'   returns = returns,
#'   W = W,
#'   m = m,
#'   bandwidth = bandwidth,
#'   kernel_func = epanechnikov_kernel,
#'   lambda = 0.1
#' )
#'
#' # Access the forecasted returns
#' head(results$forecasts)
#'
#' @export
forecast_local_factor_model <- function(
    returns,
    W,
    m,
    bandwidth,
    kernel_func = epanechnikov_kernel,
    lambda = 0.1  # Penalization parameter
) {
  T <- nrow(returns)
  p <- ncol(returns)

  # Initialize storage
  forecasts <- matrix(NA, nrow = T, ncol = p)
  est_covariances <- vector("list", T)
  residual_covariances <- vector("list", T)
  optimal_weights <- vector("list", T)

  for (tau in (W + 1):T) {
    # Compute boundary kernel weights for times 1..(tau-1)
    weight_vec <- sapply(1:(tau - 1), function(r) {
      boundary_kernel(tau, r, tau - 1, bandwidth, kernel_func)
    })
    weight_vec <- weight_vec / sum(weight_vec)  # Normalize weights
    sqrt_w <- sqrt(weight_vec)

    # Apply weights to sub-returns
    weighted_subreturns <- sweep(returns[1:(tau - 1), ], 1, sqrt_w, `*`)

    # Perform Local PCA
    pca_result <- prcomp(weighted_subreturns, center = FALSE, scale. = FALSE)
    num_factors <- min(m, ncol(pca_result$x))
    if (num_factors < 1) next  # Skip if no factors

    # Extract factors and loadings
    Fhat <- pca_result$x[, 1:num_factors, drop = FALSE] / sqrt(tau - 1)  # Normalize
    loadings_hat <- pca_result$rotation[, 1:num_factors, drop = FALSE]

    # Forecast for time tau: Use last factor scores
    F_tau_hat <- Fhat[nrow(Fhat), , drop = FALSE]
    R_hat_tau <- loadings_hat %*% t(F_tau_hat)
    forecasts[tau, ] <- as.vector(R_hat_tau)

    # Estimate Factor Covariance
    factor_cov <- cov(Fhat, use = "complete.obs")

    # Estimate Residual Covariance
    residuals_sub <- returns[1:(tau - 1), ] - Fhat %*% t(loadings_hat)
    resid_cov <- estimate_residual_cov(residuals_sub, lambda)
    residual_covariances[[tau]] <- resid_cov

    # Combine to get Estimated Covariance
    est_cov <- loadings_hat %*% factor_cov %*% t(loadings_hat) + resid_cov
    est_covariances[[tau]] <- est_cov

    # Compute Optimal Weights using QP
    weights <- compute_optimal_weights(est_cov, p)
    optimal_weights[[tau]] <- weights
  }

  return(list(
    forecasts = forecasts,
    est_covariances = est_covariances,
    residual_covariances = residual_covariances,
    optimal_weights = optimal_weights
  ))
}
#' @keywords internal
compute_optimal_weights <- function(cov_matrix, p) {
  if (any(is.na(cov_matrix)) || !is.positive.definite(round(cov_matrix, 10))) {
    # Add a small ridge term to stabilize inversion
    cov_matrix <- cov_matrix + diag(1e-6, p)
  }

  Dmat <- solve(cov_matrix)
  dvec <- rep(0, p)
  Amat <- matrix(1, nrow = p, ncol = 1)  # Constraint: sum(weights) = 1
  bvec <- 1

  # Solve QP
  result <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    return(result$solution)
  } else {
    # Return equal weights if QP fails
    return(rep(1/p, p))
  }
}
