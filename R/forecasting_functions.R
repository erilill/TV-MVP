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
    returns         ,
    initial_window  ,  # how many periods in the initial “estimation”
    rebal_period    ,  # holding window length (HT in the paper)
    max_factors     ,
    return_type    = "daily",
    kernel_func    = epanechnikov_kernel,
    rf             = NULL,
    M0             = 10, #Cov func
    rho_grid = seq(0.005, 2, length.out = 30), #Cov func
    floor_value    = 1e-12, #Cov func
    epsilon2       = 1e-6) { #Cov func
  iT <- nrow(returns)
  ip <- ncol(returns)
  rebalance_dates <- seq(initial_window + 1, iT, by = rebal_period)
  RT <- length(rebalance_dates)
  
  if (is.null(rf)) {
    rf_vec <- rep(0, (iT-initial_window))
  } else {
    rf_vec <- if (length(rf) == 1) rep(rf, (iT-initial_window)) else rf
  }
  
  # Initialize storage
  weights <- list()
  cum_rebal_returns <- numeric(RT)
  daily_port_ret <- numeric(0)
  theoretical_risk <- numeric(0) # sqrt(w' Sigma w)
  theoretical_mu  <- numeric(0) 
  
  # Determine number of factors <- would be good to re-compute yearly
  m <- determine_factors(returns[1:initial_window,], max_factors, silverman(returns[1:initial_window,]))$optimal_R
  
  for (l in seq_len(RT)) {
    reb_t <- rebalance_dates[l]
    est_data <- returns[1:(reb_t - 1), , drop=FALSE]
    bandwidth <- silverman(returns)
    
    # Local PCA
    local_res <- localPCA(est_data, bandwidth, m, kernel_func)
    
    # Compute covariance
    Sigma_hat <- estimate_residual_cov_poet_local(localPCA_results =  local_res,
                                                  returns = est_data,
                                                  M0 = 10, 
                                                  rho_grid = seq(0.005, 2, length.out = 30),
                                                  floor_value = 1e-12,
                                                  epsilon2 = 1e-6)$total_cov
    
    # Compute weights
    inv_cov <- chol2inv(chol(Sigma_hat))
    ones <- rep(1, ip)
    w_gmv_unnorm <- inv_cov %*% ones
    w_hat <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))  # Normalize weights
    
    weights[[l]] <- w_hat
    
    hold_end <- min(reb_t + rebal_period - 1, iT)
    port_ret_window <- returns[reb_t:hold_end, , drop=FALSE] %*% w_hat
    
    L <- hold_end - reb_t+1
    
    # Daily realized returns
    daily_port_ret <- c(daily_port_ret, port_ret_window)
    if (l == 1) {
      cum_rebal_returns[l] <- sum(port_ret_window)
    } else {
      cum_rebal_returns[l] <- cum_rebal_returns[l - 1] + sum(port_ret_window)
    }
  }
  
  
  # Cumulative returns
  N <- length(daily_port_ret)
  excess_ret <- daily_port_ret - rf_vec
  CER <- sum(excess_ret)
  
  # Metrics
  mean_val <- CER / N
  sample_sd <- sd(excess_ret)
  sample_SR <- mean(excess_ret) / sample_sd
  
  # Set annualization factor based on return frequency
  annualization_factor <- switch(return_type,
                                 "daily" = sqrt(252),  # Daily returns
                                 "monthly" = sqrt(12), # Monthly returns
                                 "weekly" = sqrt(52),  # Weekly returns
                                 stop("Invalid return type! Choose 'daily', 'monthly', or 'weekly'.")
  )
  
  sample_sd_annualized <- sample_sd*annualization_factor #wrong?
  sample_SR_annualized <- sample_SR*annualization_factor # wrong?
  
  
  list(
    rebal_dates              = rebalance_dates,
    weights                  = weights,
    excess_returns           = excess_ret,
    cum_rebal_returns        = cum_rebal_returns,
    cumulative_excess_return = CER,
    standard_deviation       = sample_sd,
    sharpe_ratio             = sample_SR,
    standard_deviation_annualized = sample_sd_annualized,
    sharpe_ratio_annualized = sample_SR_annualized
  )
}
#' Predict Optimal Portfolio Weights Using Time-Varying Covariance Estimation
#'
#' This function estimates the Global Minimum Variance Portfolio (MVP)
#' and optionally a Minimum Variance Portfolio (MVP) with a return constraint, 
#' using time-varying covariance estimation based on Local Principal Component Analysis (Local PCA).
#' 
#' @param returns A numeric matrix of asset returns (T x p), where T is the number of time periods and p is the number of assets.
#' @param horizon Integer. Investment horizon over which expected return and risk are computed. Default is `1`.
#' @param max_factors Integer. The maximum number of latent factors to consider in the Local PCA model. Default is `3`.
#' @param kernel_func Function. A kernel function used for weighting observations in Local PCA. Default is `epanechnikov_kernel`.
#' @param min_return Optional numeric. If provided, the function also computes a Minimum Variance Portfolio
#' that meets this minimum expected return constraint.
#' @param rf The risk-free rate (scalar). 
#' If `NULL`, it is assumed to be `0`. Default is `NULL`.
#'
#' @return A list with at least one component:
#' - **`GMV`**: A list containing:
#'     - `weights`: The optimal MVP weights (vector of length `p`).
#'     - `expected_return`: Expected return of the GMVP.
#'     - `risk`: Standard deviation (risk) of the GMVP.
#' 
#' If `min_return` is specified, an additional component is returned:
#' - **`MinVarWithReturnConstraint`**: A list containing:
#'     - `weights`: The optimal MVP weights under the return constraint.
#'     - `expected_return`: Expected return of the constrained MVP.
#'     - `risk`: Standard deviation (risk) of the constrained MVP.
#'
#' @details
#' The function implements a time-varying PCA approach to estimate latent factor structures 
#' and uses a sparse residual covariance estimation method to improve covariance matrix estimation.
#' The covariance matrix is used to determine the global minimum variance portfolio (MVP). 
#' The number of factors is determined by a BIC-type information criterion using the function 
#' `determine_factors`, and the bandwidth is determined by Silverman's rule of thumb. If 
#' `min_return` is provided, solves a quadratic optimization problem to find an MVP 
#' that satisfies the expected return constraint.
#' 
#' @examples
#' \dontrun{
#' # Simulate random asset returns (200 time periods, 20 assets)
#' set.seed(123)
#' returns <- matrix(rnorm(200 * 20, mean = 0, sd = 0.02), ncol = 20)
#'
#' # Compute GMVP
#' result <- predict_portfolio(returns)
#' print(result$GMV$weights)
#'
#' # Compute GMVP with return constraint
#' result <- predict_portfolio(returns, min_return = 0.005)
#' print(result$MinVarWithReturnConstraint$weights)
#'}
#'
#' @export
predict_portfolio <- function(
    returns,
    horizon = 1,
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    min_return = NULL,
    rf = NULL
) {
  iT <- nrow(returns)
  ip <- ncol(returns)
  
  # Determine optimal number of factors using Silverman’s bandwidth
  m <- determine_factors(returns, max_factors, silverman(returns))$optimal_R
  
  # Select bandwidth
  bandwidth <- silverman(returns)
  
  # Local PCA
  local_res <- localPCA(returns, bandwidth, m, kernel_func)
  
  # Compute covariance
  Sigma_hat <- estimate_residual_cov_poet_local(localPCA_results =  local_res,
                                                returns = returns,
                                                M0 = 10, 
                                                rho_grid = seq(0.005, 2, length.out = 30),
                                                floor_value = 1e-12,
                                                epsilon2 = 1e-6)$total_cov
  
  # Expected returns
  if (is.null(rf)) {
    mean_returns <- colMeans(returns)
  } else {
    mean_returns <- colMeans(returns) - rf
  }
  
  
  ## Global Minimum Variance Portfolio (GMVP)
  inv_cov <- chol2inv(chol(Sigma_hat))
  ones <- rep(1, ip)
  w_gmv_unnorm <- inv_cov %*% ones
  w_gmv <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))  # Normalize weights
  
  ## Compute GMVP Expected Return and Risk
  expected_return_gmv <- sum(w_gmv * mean_returns) * horizon
  risk_gmv <- sqrt(as.numeric(t(w_gmv) %*% Sigma_hat %*% w_gmv)) * sqrt(horizon)
  
  ### **Minimum Variance Portfolio with Return Constraint**
  if (!is.null(min_return)) {
    A  <- cbind(rep(1, ip), mean_returns)  # Constraints matrix (p x 2)
    b  <- c(1, min_return / horizon)  # Constraint values
    
    A_Sigma_inv_A <- solve(t(A) %*% inv_cov %*% A)
    
    w_constrained <- inv_cov %*% A %*% A_Sigma_inv_A %*% b
    
    # Compute Expected Return and Risk for Constrained Portfolio
    expected_return_constrained <- sum(w_constrained * mean_returns) * horizon
    risk_constrained <- sqrt(as.numeric(t(w_constrained) %*% Sigma_hat %*% w_constrained)) * sqrt(horizon)
    
    return(list(
      GMV = list(
        weights = w_gmv,
        expected_return = expected_return_gmv,
        risk = risk_gmv
      ),
      MinVarWithReturnConstraint = list(
        weights = w_constrained,
        expected_return = expected_return_constrained,
        risk = risk_constrained
      )
    ))
  } else {
    # If no return constraint, return only GMV
    return(list(
      GMV = list(
        weights = w_gmv,
        expected_return = expected_return_gmv,
        risk = risk_gmv
      )
    ))
  }
}
