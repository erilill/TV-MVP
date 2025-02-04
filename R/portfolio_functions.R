# R/portfolio_functions.R

#' Compute Optimal Portfolio Weights
#'
#' This function computes the optimal portfolio weights based on time-varying covariance matrices
#' using quadratic programming.
#'
#' @param W An integer specifying the window size for portfolio optimization.
#' @param factors_list A list of factor matrices for each time period.
#' @param loadings_list A list of loading matrices corresponding to each factor matrix.
#' @param residual_covariance A covariance matrix of the residuals.
#' @param T_total An integer representing the total number of time periods.
#' @param p An integer representing the number of assets.
#'
#' @return A list of optimal portfolio weights for each time period beyond the window.
#'
#' @details
#' The function iterates over each time period beyond the initial window \code{W} and computes
#' the optimal portfolio weights that minimize variance using the inverse of the time-varying
#' covariance matrix. It utilizes the \code{quadprog} package for quadratic programming.
#'
#' @examples
#' # Example parameters
#' W <- 60
#' factors_list <- list(matrix(rnorm(300), ncol = 5))
#' loadings_list <- list(matrix(runif(25), ncol = 5))
#' residual_covariance <- diag(0.01, 5)
#' T_total <- 200
#' p <- 5
#'
#' # Compute optimal weights
#' optimal_weights <- optimal_weights(W, factors_list, loadings_list, residual_covariance, T_total, p)
#' print(optimal_weights[[61]])
#'
#' @import quadprog
#' @export
optimal_weights <- function(W, factors_list, loadings_list, residual_covariance, T_total, p){
  optimal_weights <- list()
  for (t in (W + 1):T_total) {
    time_varying_cov <- loadings_list[[t]] %*% cov(factors_list[[t]]) %*% t(loadings_list[[t]]) + residual_covariance[[t]]
    Dmat <- solve(time_varying_cov)  # Inverse covariance matrix
    dvec <- rep(0, p)
    Amat <- cbind(rep(1, p))  # Constraint: sum of weights = 1
    bvec <- 1

    # Solve quadratic programming problem
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
    optimal_weights[[t]] <- result$solution
  }
  return(optimal_weights)
}
#'
#'
#' Compute Realized Covariance Matrices
#'
#' This function computes realized covariance matrices and residual covariance matrices
#' over a specified rolling window.
#'
#' @param W An integer specifying the window size for computing realized covariances.
#' @param returns A matrix of asset returns with rows representing time periods and
#' columns representing assets.
#' @param residuals A matrix of residuals corresponding to asset returns.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{\code{realized_covariances}}{A list of realized covariance matrices for each
#'   time period beyond the initial window.}
#'   \item{\code{realized_resid_covariances}}{A list of realized residual covariance
#'   matrices for each time period beyond the initial window.}
#' }
#'
#' @details
#' For each time period \eqn{t} beyond the initial window \code{W}, the function computes:
#' \enumerate{
#'   \item The realized covariance matrix using the returns from time \eqn{t-W} to \eqn{t}.
#'   \item The realized residual covariance matrix using the residuals from time \eqn{t-W} to \eqn{t}.
#' }
#'
#' If residuals contain missing values (\code{NA}), the corresponding residual covariance
#' matrix is set to \code{NA}.
#'
#' @examples
#' # Example parameters
#' W <- 60
#' T <- 200
#' p <- 5
#'
#' # Simulate returns and residuals
#' returns <- matrix(rnorm(T * p, mean = 0.01, sd = 0.02), ncol = p)
#' residuals_matrix <- matrix(rnorm(T * p, mean = 0, sd = 0.01), ncol = p)
#'
#' # Compute realized covariances
#' cov_results <- realized_covariances(W, returns, residuals_matrix)
#' print(cov_results$realized_covariances[[W + 1]])
#'
#' @export
realized_covariances <- function(W, returns, residuals){
  T <- nrow(returns)
  realized_covariances <- vector("list", T)
  realized_resid_covariances <- vector("list", T)
  for (t in (W + 1):T) {
    # Compute Realized Covariance
    realized_covariances[[t]] <- cov(returns[(t - W):t, ])

    # Compute Realized Residual Covariance
    # Ensure residuals are available and complete
    residuals_window <- residuals[(t - W):t, ]
    if (any(is.na(residuals_window))) {
      realized_resid_covariances[[t]] <- NA
    } else {
      realized_resid_covariances[[t]] <- cov(residuals_window)
    }
  }
  return(list(realized_covariances = realized_covariances,
              realized_resid_covariances = realized_resid_covariances))
}
#' Compute Realized Sharpe Ratios
#'
#' This function computes the realized Sharpe Ratios for each time period beyond
#' the initial window based on portfolio returns and risks.
#'
#' @param W An integer specifying the window size for computing Sharpe Ratios.
#' @param returns A matrix of asset returns with rows representing time periods and
#' columns representing assets.
#' @param loadings_list A list of loading matrices corresponding to each time period.
#' @param factors_list A list of factor matrices corresponding to each time period.
#' @param optimal_weights A list of optimal portfolio weights for each time period.
#' @param factor_covariance_matrix A covariance matrix of the factors.
#' @param residual_covariance_matrix A covariance matrix of the residuals.
#' @param risk_free_rate A numeric scalar specifying the risk-free rate. Defaults to \code{0}.
#'
#' @return A numeric vector of realized Sharpe Ratios for each time period beyond
#' the initial window.
#'
#' @details
#' For each time period \eqn{t} beyond the initial window \code{W}, the function computes:
#' \enumerate{
#'   \item The realized returns of the portfolio over the window \eqn{t-W} to \eqn{t}.
#'   \item The predicted portfolio standard deviation using the time-varying covariance matrix.
#'   \item The Sharpe Ratio as the ratio of mean portfolio returns to portfolio risk.
#' }
#'
#' @examples
#' # Example parameters
#' W <- 60
#' T <- 200
#' p <- 5
#'
#' # Simulate returns, factors, and loadings
#' returns <- matrix(rnorm(T * p, mean = 0.01, sd = 0.02), ncol = p)
#' residuals_matrix <- matrix(rnorm(T * p, mean = 0, sd = 0.01), ncol = p)
#' factors_list <- list(matrix(rnorm(300), ncol = 5))
#' loadings_list <- list(matrix(runif(25), ncol = 5))
#'
#' # Compute residual covariance
#' residual_covariance <- estimate_residual_cov(residuals_matrix)
#'
#' # Compute time-varying covariance
#' factor_cov <- factor_covariance(factors_list)
#' Omega_latest <- compute_time_varying_cov(loadings_list[[1]], factor_cov, residual_covariance)
#'
#' # Compute optimal weights
#' optimal_weights <- compute_optimal_weights(W, factors_list, loadings_list, residual_covariance, T, p)
#'
#' # Compute realized covariances
#' cov_results <- realized_covariances(W, returns, residuals_matrix)
#'
#' # Compute Sharpe Ratios
#' sharpe_ratios <- SR(W, returns, loadings_list, factors_list, optimal_weights)
#' print(sharpe_ratios)
#'
#' @export
SR <- function(W, returns, loadings_list, factors_list, optimal_weights, residual_covariance_matrix, risk_free_rate = 0){
  T <- nrow(returns)
  realized_sharpes <- rep(NA, T)
  for (t in (W+1):T) {
    tmp <- returns[(t - W):t, , drop = FALSE]

    # Weighted portfolio returns
    weights <- optimal_weights[[t]]
    portfolio_returns <- rowSums(sweep(tmp, 2, weights, `*`))  # Vector of portfolio returns

    # Compute mean portfolio return
    mean_rets <- mean(portfolio_returns)

    # Compute portfolio risk using time-varying covariance matrix
    time_varying_cov <- loadings_list[[t]] %*% cov(factors_list[[t]]) %*% t(loadings_list[[t]]) + residual_covariance_matrix[[t]]
    portfolio_std_dev <- sqrt(as.numeric(t(weights) %*% time_varying_cov %*% weights))

    # Compute Sharpe Ratio
    realized_sharpes[t] <- (mean_rets - risk_free_rate) / portfolio_std_dev
  }
  return(realized_sharpes)
}
#' @export
compute_mean_returns <- function(local_loadings, local_factors) {
  # Check inputs
  if (length(local_loadings) != length(local_factors)) {
    stop("The number of local loadings must match the number of local factors.")
  }

  T_total <- length(local_loadings)   # Total number of time periods
  p <- nrow(local_loadings[[1]])      # Number of assets
  m <- ncol(local_loadings[[1]])      # Number of factors

  # Ensure all loadings and factors have consistent dimensions
  for (t in 1:T_total) {
    if (!all(dim(local_loadings[[t]]) == c(p, m))) {
      stop(sprintf("Dimension mismatch in local_loadings at time %d.", t))
    }
    if (!all(dim(local_factors[[t]]) == c(1, m))) {
      stop(sprintf("Dimension mismatch in local_factors at time %d.", t))
    }
  }

  # 1. Compute the average factor vector across all time periods
  #    Resulting in a 1 x m vector
  avg_factors <- colMeans(do.call(rbind, local_factors))  # 1 x m

  # 2. Initialize a matrix to store mean returns for each time t
  #    Resulting in a T_total x p matrix
  mu_matrix <- matrix(0, nrow = T_total, ncol = p)
  colnames(mu_matrix) <- paste0("Asset_", 1:p)
  rownames(mu_matrix) <- paste0("Time_", 1:T_total)

  # 3. Compute mu_t for each time period t
  for (t in 1:T_total) {
    B_t <- local_loadings[[t]]         # p x m matrix
    mu_t <- B_t %*% avg_factors         # p x 1 vector
    mu_matrix[t, ] <- as.vector(mu_t)    # Store as a row in mu_matrix
  }

  return(mu_matrix)  # T_total x p matrix of predicted mean returns
}

