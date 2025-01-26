# R/portfolio_functions.R

#' Compute Minimum-Variance Portfolios
#'
#' This function computes the global minimum-variance portfolio and, optionally,
#' the optimal portfolio achieving a specified target return (\code{mu_p}).
#' It returns both variance and standard deviation (risk) for each portfolio.
#'
#' @param mu A numeric vector of expected returns for assets. Typically computed as
#' \code{colMeans(returns, na.rm = TRUE)}.
#' @param Omega A positive definite covariance matrix of asset returns.
#' @param mu_p An optional numeric scalar specifying the target expected return for the
#' portfolio. Defaults to \code{NULL}, indicating that only the global minimum-variance
#' portfolio is computed.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{w_g}}{Numeric vector. Weights of the global minimum-variance portfolio.}
#'   \item{\code{mu_g}}{Numeric scalar. Expected return of the global minimum-variance portfolio.}
#'   \item{\code{var_g}}{Numeric scalar. Variance of the global minimum-variance portfolio.}
#'   \item{\code{risk_g}}{Numeric scalar. Standard deviation (risk) of the global minimum-variance portfolio.}
#'   \item{\code{w_p}}{Numeric vector. Weights of the portfolio achieving the target expected return \code{mu_p}.}
#'   \item{\code{mu_p}}{Numeric scalar. Target expected return of the portfolio.}
#'   \item{\code{var_p}}{Numeric scalar. Variance of the portfolio achieving \code{mu_p}.}
#'   \item{\code{risk_p}}{Numeric scalar. Standard deviation (risk) of the portfolio achieving \code{mu_p}.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates the input covariance matrix for positive definiteness.
#'   \item Computes the global minimum-variance portfolio.
#'   \item If a target return \code{mu_p} is provided, computes the optimal portfolio
#'   achieving that return.
#' }
#'
#' The function ensures that the covariance matrix is positive definite and handles
#' any computational issues by providing informative warnings and default equal weights
#' when necessary.
#'
#' @examples
#' # Simulate data for 5 assets over 200 time periods
#' set.seed(123)
#' p <- 5  # Number of assets
#' T <- 200  # Number of time periods
#'
#' # Simulate returns matrix
#' returns <- matrix(rnorm(T * p, mean = 0.01, sd = 0.02), ncol = p)
#'
#' # Compute expected returns as column means
#' expected_returns <- colMeans(returns, na.rm = TRUE)
#'
#' # Compute covariance matrix for the latest time period (using entire returns matrix)
#' Omega_latest <- cov(returns, use = "complete.obs")
#'
#' # Compute global minimum-variance portfolio
#' global_portfolio <- minvar_portfolio(mu = expected_returns, Omega = Omega_latest)
#' print(global_portfolio)
#'
#' # Compute both global and targeted minimum-variance portfolios with mu_p = 0.05
#' targeted_portfolio <- minvar_portfolio(mu = expected_returns, Omega = Omega_latest, mu_p = 0.05)
#' print(targeted_portfolio)
#'
#' @export
minvar_portfolio <- function(mu, Omega, mu_p = NULL) {
  # Validate inputs
  if (!is.numeric(mu) || !is.vector(mu)) {
    stop("mu must be a numeric vector of expected returns.")
  }

  if (!is.matrix(Omega) || nrow(Omega) != length(mu) || ncol(Omega) != length(mu)) {
    stop("Omega must be a square covariance matrix with dimensions matching the length of mu.")
  }

  # Check if covariance matrix is positive definite
  if (!requireNamespace("matrixcalc", quietly = TRUE)) {
    stop("Package 'matrixcalc' is required but not installed. Please install it using install.packages('matrixcalc').")
  }

  if (!matrixcalc::is.positive.definite(round(Omega, 10))) {
    warning("Covariance matrix is not positive definite. Adding a small ridge term to make it positive definite.")
    Omega <- Omega + diag(1e-6, nrow = Omega)

    if (!matrixcalc::is.positive.definite(round(Omega, 10))) {
      warning("Covariance matrix remains non-positive definite after ridge adjustment. Assigning equal weights.")
      Omega <- diag(1 / length(mu), nrow = Omega)
    }
  }

  # Perform Cholesky decomposition safely
  chol_Omega <- tryCatch(chol(Omega), error = function(e) NULL)

  if (is.null(chol_Omega)) {
    warning("Cholesky decomposition failed. Assigning equal weights for global portfolio.")
    w_g <- rep(1 / length(mu), length(mu))
    mu_g <- sum(w_g * mu)
    var_g <- sum(w_g^2 * diag(Omega))
    risk_g <- sqrt(var_g)

    portfolio <- list(
      w_g = w_g,
      mu_g = mu_g,
      var_g = var_g,
      risk_g = risk_g
    )

    if (!is.null(mu_p)) {
      portfolio$w_p <- rep(1 / length(mu), length(mu))
      portfolio$mu_p <- sum(portfolio$w_p * mu)
      portfolio$var_p <- sum(portfolio$w_p^2 * diag(Omega))
      portfolio$risk_p <- sqrt(portfolio$var_p)
    }

    return(portfolio)
  }

  # Compute inverse of Omega
  inv_Omega <- chol2inv(chol_Omega)

  # Compute necessary sums
  tmpa <- as.numeric(t(mu) %*% inv_Omega %*% mu)
  tmpb <- sum(inv_Omega %*% mu)
  tmpc <- sum(inv_Omega)

  # Compute global minimum-variance portfolio
  mu_g <- tmpb / tmpc
  w_g <- as.numeric(inv_Omega %*% rep(1, length(mu)) / tmpc)
  var_g <- as.numeric(t(w_g) %*% Omega %*% w_g)
  risk_g <- sqrt(var_g)

  # Initialize portfolio list
  portfolio <- list(
    w_g = w_g,
    mu_g = mu_g,
    var_g = var_g,
    risk_g = risk_g
  )

  # If mu_p is provided, compute the optimal portfolio achieving mu_p
  if (!is.null(mu_p)) {
    # Check if mu_p is achievable
    if (mu_p < min(mu) || mu_p > max(mu)) {
      warning(sprintf("Target return mu_p = %.4f is outside the range of asset returns. Skipping optimal portfolio computation.", mu_p))
    } else {
      # Compute optimal portfolio weights
      numerator <- inv_Omega %*% (mu * (tmpc * mu_p - tmpb) + (tmpa - tmpb * mu_p))
      denominator <- tmpa * tmpc - tmpb^2

      if (denominator == 0) {
        warning("Denominator for optimal portfolio computation is zero. Assigning equal weights.")
        w_p <- rep(1 / length(mu), length(mu))
      } else {
        w_p <- as.numeric(numerator / denominator)
      }

      # Compute variance and risk of the optimal portfolio
      var_p <- as.numeric(t(w_p) %*% Omega %*% w_p)
      risk_p <- sqrt(var_p)

      # Add to portfolio list
      portfolio$w_p <- w_p
      portfolio$mu_p <- mu_p
      portfolio$var_p <- var_p
      portfolio$risk_p <- risk_p
    }
  }

  # Return the portfolio list
  return(portfolio)
}
#' Plot Portfolio Weights
#'
#' This function creates a bar plot of portfolio weights for visual analysis.
#'
#' @param weights A numeric vector representing the weights of assets in the portfolio.
#' @param title A character string specifying the title of the plot. Defaults to "Portfolio Weights".
#'
#' @return A ggplot object displaying the portfolio weights as a bar chart.
#'
#' @details
#' The function takes the asset weights and generates a bar plot using \code{ggplot2}.
#' Assets are ordered by their weights in descending order for better visualization.
#'
#' @examples
#' # Example weights
#' weights <- c(0.2, 0.3, 0.25, 0.15, 0.1)
#'
#' # Plot portfolio weights
#' plot_portfolio_weights(weights, title = "Sample Portfolio Weights")
#'
#' @import ggplot2
#' @export
plot_portfolio_weights <- function(weights, title = "Portfolio Weights") {
  df <- data.frame(
    Asset = paste0("Asset ", 1:length(weights)),
    Weight = weights
  )

  ggplot(df, aes(x = reorder(Asset, -Weight), y = Weight)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = title, x = "Asset", y = "Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
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
