#' @export
rolling_time_varying_mvp <- function(
    returns         ,
    initial_window  ,  # how many periods in the initial “estimation”
    rebal_period    ,  # holding window length (HT in the paper)
    max_factors     ,
    return_type    = "daily",
    kernel_func    = epanechnikov_kernel,
    rf             = NULL ) {
  T <- nrow(returns)
  p <- ncol(returns)
  rebalance_dates <- seq(initial_window + 1, T, by = rebal_period)
  RT <- length(rebalance_dates)
  
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
    Sigma_hat <- estimate_residual_cov_poet_local(local_res, est_data)$total_cov
    
    # Compute weights
    inv_cov <- chol2inv(chol(Sigma_hat))
    ones <- rep(1, p)
    w_gmv_unnorm <- inv_cov %*% ones
    w_hat <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))  # Normalize weights
    
    weights[[l]] <- w_hat
    
    hold_end <- min(reb_t + rebal_period - 1, T)
    port_ret_window <- returns[reb_t:hold_end, , drop=FALSE] %*% w_hat
    
    ## Theoretical metrics
    # Expected returns
    mean_returns <- colMeans(est_data) # place holder, might change
    
    ## Compute Expected Return and Risk
    mu <- sum(w_hat * mean_returns)
    risk <- sqrt(as.numeric(t(w_hat) %*% Sigma_hat %*% w_hat))
    L <- hold_end - reb_t+1
    theoretical_mu <- c(theoretical_mu, mu*(1:L))
    theoretical_risk <- c(theoretical_risk, risk*sqrt(1:L))
    
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
  excess_ret <- daily_port_ret - rf
  theoretical_mu <- theoretical_mu - rf
  CER <- sum(excess_ret)
  
  # Metrics
  mean_val <- CER / N
  sample_sd <- sd(excess_ret)
  sample_SR <- mean(excess_ret) / sample_sd
  RMSE <- sqrt(mean((excess_ret - (theoretical_mu))^2))  
  avg_risk_diff <- abs(sample_sd - mean(theoretical_risk))
  
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
    RMSE                     = RMSE,
    standard_deviation_annualized = sample_sd_annualized,
    sharpe_ratio_annualized = sample_SR_annualized,
    theoretical_mu = theoretical_mu,
    theoretical_risk = theoretical_risk,
    avg.risk_diff = avg_risk_diff
  )
}
#' @export
predict_portfolio <- function(
    returns,
    horizon = 1,
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    min_return = NULL,
    rf = NULL
) {
  T <- nrow(returns)
  p <- ncol(returns)
  
  # Determine optimal number of factors using Silverman’s bandwidth
  m <- determine_factors(returns, max_factors, silverman(returns))$optimal_R
  
  # Select bandwidth
  bandwidth <- silverman(returns)
  
  # Local PCA
  local_res <- localPCA(returns, bandwidth, m, kernel_func)
  
  # Compute covariance
  Sigma_hat <- estimate_residual_cov_poet_local(local_res, returns)$total_cov
  
  # Expected returns
  mean_returns <- colMeans(returns) - rf[T] # place holder, might change
  
  ## Global Minimum Variance Portfolio (GMVP)
  inv_cov <- chol2inv(chol(Sigma_hat))
  ones <- rep(1, p)
  w_gmv_unnorm <- inv_cov %*% ones
  w_gmv <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))  # Normalize weights
  
  ## Compute GMVP Expected Return and Risk
  expected_return_gmv <- sum(w_gmv * mean_returns) * horizon
  risk_gmv <- sqrt(as.numeric(t(w_gmv) %*% Sigma_hat %*% w_gmv)) * sqrt(horizon)
  
  ### **Minimum Variance Portfolio with Return Constraint**
  if (!is.null(min_return)) {
    A  <- cbind(rep(1, p), mean_returns)  # Constraints matrix (p x 2)
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
