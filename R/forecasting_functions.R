#' @export
rolling_time_varying_mvp <- function(
    returns         ,
    initial_window  ,  # how many periods in the initial “estimation”
    rebal_period    ,  # holding window length (HT in the paper)
    max_factors     ,
    return_type    = "daily",
    kernel_func    = epanechnikov_kernel,
    bandwidth_func = silverman) {
  T <- nrow(returns)
  p <- ncol(returns)
  rebalance_dates <- seq(initial_window + 1, T, by = rebal_period)
  RT <- length(rebalance_dates)
  weights <- matrix(NA, nrow = RT, ncol = p)

  cum_rebal_returns <- numeric(RT)
  daily_port_ret <- numeric(0)
  for (l in seq_len(RT)) {
    reb_t <- rebalance_dates[l]
    est_data <- returns[1:(reb_t - 1), , drop=FALSE]

    m <- determine_factors(est_data, max_factors, silverman(est_data))$optimal_R


    if (identical(bandwidth_func, silverman)) {
      bandwidth <- silverman(est_data)
    } else {
      bandwidth <- handle_cv_bandwidth(est_data, m, seq(0.05, 0.95, 0.05), kernel_func)
    }
    

    # Local PCA
    local_res <- local_pca(est_data, nrow(est_data), bandwidth, m, kernel_func)

    # Compute covariance
    Sigma_hat <- estimate_residual_cov_poet_local(local_res, est_data)$total_cov

    # Compute weights
    inv_cov <- solve(Sigma_hat)
    ones <- rep(1, p)
    w_gmv_unnorm <- inv_cov %*% ones
    w_hat <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))  # Normalize weights

    weights[l, ] <- w_hat

    hold_end <- min(reb_t + rebal_period - 1, T)
    port_ret_window <- returns[reb_t:hold_end, , drop=FALSE] %*% w_hat

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
  CER <- sum(daily_port_ret)

  # Metrics
  mean_val <- CER / N
  devs <- daily_port_ret - mean_val
  stdev <- sqrt( sum(devs^2) / (N - 1) )
  SR <- mean(daily_port_ret) / stdev
  
  # Set annualization factor based on return frequency
  annualization_factor <- switch(return_type,
                                 "daily" = sqrt(252),  # Daily returns
                                 "monthly" = sqrt(12), # Monthly returns
                                 "weekly" = sqrt(52),  # Weekly returns
                                 stop("Invalid return type! Choose 'daily', 'monthly', or 'weekly'.")
  )
  
  stdev_annualized <- stdev*annualization_factor
  SR_annualized <- SR*annualization_factor
  

  list(
    rebal_dates              = rebalance_dates,
    weights                  = weights,
    daily_portfolio_returns  = daily_port_ret,
    cum_rebal_returns        = cum_rebal_returns,
    cumulative_excess_return = CER,
    standard_deviation       = stdev,
    sharpe_ratio             = SR,
    standard_deviation_annualized = stdev_annualized,
    sharpe_ratio_annualized = SR_annualized
  )
}
#' @export
predict_portfolio <- function(
    returns,
    horizon = 1,
    bandwidth_func = cv_bandwidth,
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    min_return = NULL
) {
  T <- nrow(returns)
  p <- ncol(returns)

  # Determine optimal number of factors using Silverman’s bandwidth
  m <- determine_factors(returns, max_factors, silverman(returns))$optimal_R

  # Select bandwidth
  if (identical(bandwidth_func, silverman)) {
    bandwidth <- silverman(returns)
  } else {
    bandwidth <- cv_bandwidth(returns, m, seq(0.05, 0.95, 0.05), kernel_func)$optimal_h
    }

  # Local PCA
  local_res <- local_pca(returns, nrow(returns), bandwidth, m, kernel_func)
  
  # Compute covariance
  Sigma_hat <- estimate_residual_cov_poet_local(local_res, returns)$total_cov

  # Expected returns
  mean_returns <- tcrossprod(local_res$loadings, local_res$f_hat) #Correct?

  ## Global Minimum Variance Portfolio (GMVP)
  inv_cov <- solve(Sigma_hat)
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
