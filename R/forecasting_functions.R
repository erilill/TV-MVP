#' #' Rolling Time-Varying Minimum Variance Portfolio Optimization
#'
#' This function performs rolling time-varying minimum variance portfolio (TV-MVP) optimization using 
#' time-varying covariance estimation based on Local Principal Component Analysis (Local PCA). The 
#' optimization is performed over a rolling window, with periodic rebalancing.
#'
#' @param returns A matrix of asset returns (T x p), where T is the number of time periods and p is the number of assets.
#' @param initial_window An integer specifying the number of periods used in the initial estimation window.
#' @param rebal_period An integer specifying the number of periods between portfolio rebalancing (e.g., weekly, monthly).
#' @param max_factors An integer indicating the maximum number of latent factors to consider in the factor model.
#' @param return_type A string indicating the return frequency. Options: `"daily"`, `"weekly"`, or `"monthly"`. Used for annualization of evaluation metrics.
#' @param kernel_func A kernel function to be used in the local PCA procedure. Default is `epanechnikov_kernel`.
#' @param rf A numeric vector or single value representing the risk-free rate. If `NULL`, it defaults to zero.
#' @param M0 An integer specifying the number of observations to leave out between the two sub-samples for the adaptive thresholding of the residual covariance estimation.
#' @param rho_grid A numeric sequence specifying the grid of rho values for threshold selection in covariance shrinkage. Default is `seq(0.005, 2, length.out = 30)`.
#' @param floor_value A small positive value to ensure numerical stability in the covariance matrix. Default is `1e-12`.
#' @param epsilon2 A small positive value used in the adaptive thresholding of the residual covariance estimation. Default is `1e-6`.
#'
#' @return A list containing:
#' \item{rebal_dates}{Vector of rebalancing dates.}
#' \item{weights}{A list of portfolio weights at each rebalancing date.}
#' \item{excess_returns}{A numeric vector of excess portfolio returns.}
#' \item{cum_rebal_returns}{A numeric vector of cumulative portfolio returns.}
#' \item{cumulative_excess_return}{Total cumulative excess return over the period.}
#' \item{standard_deviation}{Standard deviation of the excess returns.}
#' \item{sharpe_ratio}{Sharpe ratio of the strategy.}
#' \item{standard_deviation_annualized}{Annualized standard deviation of returns.}
#' \item{sharpe_ratio_annualized}{Annualized Sharpe ratio.}
#'
#' @details
#' The function implements a rolling time-varying PCA approach to estimate latent factor structures 
#' and uses a sparse residual covariance estimation method to improve covariance matrix estimation.
#' The covariance matrix is used to determine the global minimum variance portfolio (MVP), which is 
#' rebalanced periodically according to the specified `rebal_period`. The number of factors is
#' determined by a BIC-type information criterion using the function `determine_factors`, and the
#' bandwidth is determined by Silverman's rule of thumb.
#'
#' If `rf` is `NULL`, the risk-free rate is assumed to be zero.
#'
#' @examples
#' \dontrun{
#' # Generate random returns for 20 assets over 500 periods
#' set.seed(123)
#' returns <- matrix(rnorm(5000), nrow = 500, ncol = 20)
#'
#' # Run rolling TV-MVP optimization
#' results <- rolling_time_varying_mvp(
#'   returns = returns,
#'   initial_window = 100,
#'   rebal_period = 20,
#'   max_factors = 3,
#'   return_type = "daily",
#'   kernel_func = epanechnikov_kernel,
#'   rf = NULL
#' )
#'
#' # View Sharpe ratio
#' print(results$sharpe_ratio)
#' }
#'
#' @export
rolling_time_varying_mvp <- function(
    returns,            # Log returns matrix
    initial_window,     # Number of periods in the initial estimation window
    rebal_period,       # Holding window length (HT in the paper)
    max_factors,
    return_type    = "daily",
    kernel_func    = epanechnikov_kernel,
    rf             = NULL,
    M0             = 10,                      # For covariance function
    rho_grid       = seq(0.005, 2, length.out = 30),  # For covariance function
    floor_value    = 1e-12,                   # For covariance function
    epsilon2       = 1e-6                     # For covariance function
) {
  iT <- nrow(returns)
  ip <- ncol(returns)
  rebalance_dates <- seq(initial_window + 1, iT, by = rebal_period)
  RT <- length(rebalance_dates)
  
  # Set risk-free rate vector: length should match the number of return observations after initial_window.
  if (is.null(rf)) {
    rf_vec <- rep(0, (iT - initial_window))
  } else {
    rf_vec <- if (length(rf) == 1) rep(rf, (iT - initial_window)) else rf
  }
  
  # Initialize storage for portfolio returns and weights
  tvmvp_weights <- list()
  daily_port_ret_tvmvp <- numeric(0)
  daily_port_ret_equal <- numeric(0)
  
  # Determine number of factors using the initial window
  m <- determine_factors(returns[1:initial_window, ], max_factors, 
                         silverman(returns[1:initial_window, ]))$optimal_m
  
  for (l in seq_len(RT)) {
    reb_t <- rebalance_dates[l]
    est_data <- returns[1:(reb_t - 1), , drop = FALSE]
    # Use the estimation window to compute the bandwidth
    bandwidth <- silverman(est_data)
    
    ## TVMVP Portfolio: Local PCA and Covariance Estimation
    local_res <- localPCA(est_data, bandwidth, m, kernel_func)
    Sigma_hat <- estimate_residual_cov_poet_local(
      localPCA_results = local_res,
      returns = est_data,
      M0 = M0, 
      rho_grid = rho_grid,
      floor_value = floor_value,
      epsilon2 = epsilon2
    )$total_cov
    
    # Compute weights for the minimum variance portfolio (TVMVP)
    inv_cov <- chol2inv(chol(Sigma_hat))
    ones <- rep(1, ip)
    w_gmv_unnorm <- inv_cov %*% ones
    w_hat <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))
    tvmvp_weights[[l]] <- w_hat
    
    # Determine the holding period
    hold_end <- min(reb_t + rebal_period - 1, iT)
    
    # Portfolio returns over the holding period:
    # TVMVP Portfolio returns
    port_ret_window_tvmvp <- as.numeric(returns[reb_t:hold_end, , drop = FALSE] %*% w_hat)
    daily_port_ret_tvmvp <- c(daily_port_ret_tvmvp, port_ret_window_tvmvp)
    
    # Equal Weights Portfolio returns (baseline)
    w_equal <- rep(1/ip, ip)
    port_ret_window_equal <- as.numeric(returns[reb_t:hold_end, , drop = FALSE] %*% w_equal)
    daily_port_ret_equal <- c(daily_port_ret_equal, port_ret_window_equal)
  }
  
  ## Compute performance metrics for both portfolios
  
  # TVMVP Portfolio
  excess_ret_tvmvp <- daily_port_ret_tvmvp - rf_vec
  CER_tvmvp <- sum(excess_ret_tvmvp)
  sd_tvmvp <- sqrt(var(excess_ret_tvmvp))
  mean_ret_tvmvp <- mean(excess_ret_tvmvp)
  SR_tvmvp <- mean_ret_tvmvp / sd_tvmvp
  
  # Equal Weights Portfolio
  excess_ret_equal <- daily_port_ret_equal - rf_vec
  CER_equal <- sum(excess_ret_equal)
  sd_equal <- sqrt(var(excess_ret_equal))
  mean_ret_equal <- mean(excess_ret_equal)
  SR_equal <- mean_ret_equal / sd_equal
  
  # Set annualization factor based on return_type
  annualization_factor <- switch(return_type,
                                 "daily" = sqrt(252),
                                 "monthly" = sqrt(12),
                                 "weekly" = sqrt(52),
                                 stop("Invalid return type! Choose 'daily', 'monthly', or 'weekly'.")
  )
  
  mean_annualized_tvmvp <- mean_ret_tvmvp * (annualization_factor^2)
  sd_annualized_tvmvp   <- sd_tvmvp * annualization_factor
  
  mean_annualized_equal <- mean_ret_equal * (annualization_factor^2)
  sd_annualized_equal   <- sd_equal * annualization_factor
  
  # Compile a summary table of results
  summary_df <- data.frame(
    Method = c("Time-Varying MVP", "Equal Weight"),
    Cumulative_Excess_Return = c(CER_tvmvp, CER_equal),
    Mean_Excess_Return       = c(mean_ret_tvmvp, mean_ret_equal),
    Standard_Deviation       = c(sd_tvmvp, sd_equal),
    Sharpe_Ratio             = c(SR_tvmvp, SR_equal),
    Mean_Annualized          = c(mean_annualized_tvmvp, mean_annualized_equal),
    SD_Annualized            = c(sd_annualized_tvmvp, sd_annualized_equal)
  )
  
  # -- Construct named lists for the TVMVP and Equal strategies:
  TVMVP <- list(
    rebal_dates = rebalance_dates,
    weights     = tvmvp_weights,
    returns     = excess_ret_tvmvp
  )
  
  Equal <- list(
    rebal_dates = rebalance_dates,
    weights     = rep(1 / ip, ip),  # or you can store them by rebalancing period if needed
    returns     = excess_ret_equal
  )
  
  # Create and return an R6 object:
  out <- RollingWindow$new(
    summary = summary_df,
    TVMVP   = TVMVP,
    Equal   = Equal
  )
  
  return(out)
}
#' Predict Optimal Portfolio Weights Using Time-Varying Covariance Estimation
#'
#' This function estimates optimal portfolio weights by employing a time-varying covariance estimation
#' approach based on Local Principal Component Analysis (Local PCA). It computes the following portfolios:
#' \enumerate{
#'   \item The **Minimum Variance Portfolio (GMV)**, i.e., the portfolio with the lowest variance.
#'   \item The **Maximum SR Portfolio**, i.e., the portfolio that maximizes the Sharpe ratio (if \code{max_SR = TRUE}).
#'   \item The **Return-Constrained Portfolio**, i.e., the portfolio that minimizes variance while meeting a specified minimum expected return (if \code{min_return} is provided).
#' }
#'
#' @param returns A numeric matrix of asset returns (T x p), where T is the number of time periods and p is the number of assets.
#' @param horizon Integer. Investment horizon over which expected return and risk are computed. Default is 1.
#' @param max_factors Integer. The maximum number of latent factors to consider in the Local PCA model. Default is 3.
#' @param kernel_func Function. A kernel function used for weighting observations in Local PCA. Default is \code{epanechnikov_kernel}.
#' @param min_return Optional numeric. If provided, the function computes a Return-Constrained Portfolio that meets this minimum expected return constraint.
#' @param max_SR Optional logical. If set to TRUE, the function also computes the Maximum SR Portfolio. Default is \code{NULL}.
#' @param rf Numeric scalar. The risk-free rate. If \code{NULL}, it is assumed to be 0. Default is \code{NULL}.
#'
#' @return An object of class \code{PortfolioPredictions} containing:
#' \describe{
#'   \item{summary}{A data frame summarizing the key performance metrics for each computed portfolio. The columns include:
#'      \describe{
#'         \item{Method}{Portfolio type: "Minimum Variance Portfolio", "Maximum SR Portfolio", or "Return-Constrained Portfolio".}
#'         \item{expected_return}{The portfolio's expected return.}
#'         \item{risk}{The portfolio's risk (standard deviation).}
#'         \item{sharpe}{The portfolio's Sharpe ratio.}
#'      }
#'   }
#'   \item{GMV}{A list containing the weights, expected return, risk, and Sharpe ratio for the Minimum Variance Portfolio.}
#'   \item{max_SR}{A list containing the corresponding metrics for the Maximum SR Portfolio (if computed).}
#'   \item{MinVarWithReturnConstraint}{A list containing the corresponding metrics for the Return-Constrained Portfolio (if computed).}
#' }
#'
#' @details
#' The function employs a Local PCA approach to estimate latent factor structures and uses this to obtain a time-varying
#' estimate of the covariance matrix via POET. The number of factors is determined by an information criterion (using
#' \code{determine_factors}) and the bandwidth is set according to Silverman's rule of thumb.
#'
#' Expected returns for each asset are forecasted via a simple ARIMA model selection procedure. If a risk-free rate (\code{rf})
#' is provided, it is subtracted from the forecasts to yield excess returns.
#'
#' The optimization then proceeds as follows:
#' \enumerate{
#'   \item The Global Minimum Variance Portfolio (GMV) is computed in closed form.
#'   \item If \code{max_SR = TRUE}, the Maximum Sharpe Ratio Portfolio is computed using the closed-form
#'         solution \eqn{w \propto \Sigma^{-1} (\mu - r_f)}.
#'   \item If \code{min_return} is provided, a quadratic optimization problem is solved to obtain a Return-Constrained Portfolio.
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate random asset returns (200 time periods, 20 assets)
#' set.seed(123)
#' returns <- matrix(rnorm(200 * 20, mean = 0, sd = 0.02), ncol = 20)
#'
#' # Compute portfolios using the function
#' result <- predict_portfolio(returns, horizon = 5, max_factors = 3,
#'                              min_return = 0.005, max_SR = TRUE)
#'
#' # Print the summary of portfolio performance
#' print(result)
#'
#' # Access the weights of the Minimum Variance Portfolio
#' result$GMV$weights
#' }
#'
#' @export

predict_portfolio <- function(
    returns,
    horizon = 1,
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    min_return = NULL,
    max_SR = NULL,  # flag: if TRUE, compute maximum Sharpe portfolio
    rf = NULL
) {
  iT <- nrow(returns)
  ip <- ncol(returns)
  
  # Determine optimal number of factors using Silverman’s bandwidth
  m <- determine_factors(returns, max_factors, silverman(returns))$optimal_m
  
  # Select bandwidth
  bandwidth <- silverman(returns)
  
  # Local PCA
  local_res <- localPCA(returns, bandwidth, m, kernel_func)
  
  # Compute covariance matrix from local PCA results using POET
  Sigma_hat <- estimate_residual_cov_poet_local(localPCA_results = local_res,
                                                returns = returns,
                                                M0 = 10, 
                                                rho_grid = seq(0.005, 2, length.out = 30),
                                                floor_value = 1e-12,
                                                epsilon2 = 1e-6)$total_cov
  

  # Expected returns
  if (is.null(rf)) {
    mean_returns <- comp_expected_returns(returns, horizon)
  } else {
    mean_returns <- comp_expected_returns(returns, horizon) - rf
  }
  
  ## Compute GMVP
  inv_cov <- chol2inv(chol(Sigma_hat))
  ones <- rep(1, ip)
  w_gmv_unnorm <- inv_cov %*% ones
  w_gmv <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))
  
  expected_return_gmv <- sum(w_gmv * mean_returns) * horizon
  risk_gmv <- sqrt(as.numeric(t(w_gmv) %*% Sigma_hat %*% w_gmv)) * sqrt(horizon)
  
  GMV <- list(
    weights = w_gmv,
    expected_return = expected_return_gmv,
    risk = risk_gmv,
    sharpe = (expected_return_gmv/horizon) / (risk_gmv/sqrt(horizon))
  )
  
  # Compute Maximum Sharpe Ratio portfolio if requested
  max_sr_portfolio <- NULL
  if (!is.null(max_SR) && max_SR == TRUE) {
    w_unnorm <- inv_cov %*% (mean_returns / horizon)
    w_max_sr <- as.numeric(w_unnorm / sum(w_unnorm))
    expected_return_max_sr <- sum(w_max_sr * mean_returns) * horizon
    risk_max_sr <- sqrt(as.numeric(t(w_max_sr) %*% Sigma_hat %*% w_max_sr)) * sqrt(horizon)
    
    max_sr_portfolio <- list(
      weights = w_max_sr,
      expected_return = expected_return_max_sr,
      risk = risk_max_sr,
      sharpe = (expected_return_max_sr/horizon) / (risk_max_sr/sqrt(horizon))
    )
  }
  
  # Compute Minimum Variance Portfolio with Return Constraint (if applicable)
  constrained_portfolio <- NULL
  if (!is.null(min_return)) {
    A  <- cbind(rep(1, ip), mean_returns)  # (p x 2) constraints
    b  <- c(1, min_return / horizon)
    A_Sigma_inv_A <- solve(t(A) %*% inv_cov %*% A)
    w_constrained <- inv_cov %*% A %*% A_Sigma_inv_A %*% b
    expected_return_constrained <- sum(w_constrained * mean_returns) * horizon
    risk_constrained <- sqrt(as.numeric(t(w_constrained) %*% Sigma_hat %*% w_constrained)) * sqrt(horizon)
    
    constrained_portfolio <- list(
      weights = w_constrained,
      expected_return = expected_return_constrained,
      risk = risk_constrained,
      sharpe = (expected_return_constrained/horizon) / (risk_constrained/sqrt(horizon))
    )
  }
  
  # Build summary data frame
  method_names <- c("GMV")
  expected_returns_vec <- c(expected_return_gmv)
  risk_vec <- c(risk_gmv)
  sharpe_vec <- c(GMV$sharpe)
  
  if (!is.null(max_sr_portfolio)) {
    method_names <- c(method_names, "max_SR")
    expected_returns_vec <- c(expected_returns_vec, expected_return_max_sr)
    risk_vec <- c(risk_vec, risk_max_sr)
    sharpe_vec <- c(sharpe_vec, max_sr_portfolio$sharpe)
  }
  if (!is.null(constrained_portfolio)) {
    method_names <- c(method_names, "MinVarWithReturnConstraint")
    expected_returns_vec <- c(expected_returns_vec, expected_return_constrained)
    risk_vec <- c(risk_vec, risk_constrained)
    sharpe_vec <- c(sharpe_vec, constrained_portfolio$sharpe)
  }
  
  summary_df <- data.frame(
    Method = method_names,
    expected_return = expected_returns_vec,
    risk = risk_vec,
    sharpe = sharpe_vec
  )
  
  # Create and return an object of PortfolioPredictions (an R6 object)
  out <- PortfolioPredictions$new(
    summary = summary_df,
    GMV = GMV,
    max_SR = if (!is.null(max_sr_portfolio)) max_sr_portfolio else NULL,
    MinVarWithReturnConstraint = if (!is.null(constrained_portfolio)) constrained_portfolio else NULL
  )
  
  return(out)
}

#' Function to compute expected returns using a simple model selection approach
#' @param returns T times p matrix of returns
#' @param horizon Length of forecasting horizon
comp_expected_returns <- function(returns, horizon) {
  exp_ret <- numeric(ncol(returns))
  
  for (i in seq_len(ncol(returns))) {
    candidate_models <- list()
    aics <- numeric()
    
    for (order in list(c(0,0,0), c(1,0,0), c(0,0,1), c(1,0,1))) {
      model <- tryCatch(
        arima(returns[, i], order = order),
        error = function(e) NULL
      )
      candidate_models <- c(candidate_models, list(model))
      aics <- c(aics, if (!is.null(model)) AIC(model) else Inf)
    }
    
    # If all models failed, fallback to mean return
    if (all(is.infinite(aics))) {
      exp_ret[i] <- mean(returns[, i])
    } else {
      best_model <- candidate_models[[which.min(aics)]]
      fc <- predict(best_model, n.ahead = horizon)$pred
      exp_ret[i] <- mean(fc)
    }
  }
  
  return(exp_ret)
}
