#' @export
rolling_time_varying_mvp <- function(
    returns,
    initial_window = 60,  # how many periods in the initial “estimation”
    rebal_period   = 20,  # holding window length (HT in the paper)
    max_factors    = 3,
    bandwidth      = 0.2,
    kernel_func    = epanechnikov_kernel
) {
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

    local_res <- localPCA(est_data, bandwidth, max_factors, kernel_func)
    factor_cov   <- cov(local_res$factors)
    residuals    <- compute_residuals(local_res$factors, local_res$loadings, est_data)
    residual_cov <- estimate_residual_cov(residuals)

    last_t <- nrow(est_data)
    loadings_mid <- local_res$loadings[[last_t]]

    Sigma_hat <- loadings_mid %*% factor_cov %*% t(loadings_mid) + residual_cov

    inv_cov <- solve(Sigma_hat)
    ones <- rep(1, p)
    w_unnorm <- inv_cov %*% ones
    w_hat <- as.numeric(w_unnorm / sum(w_unnorm))

    weights[l, ] <- w_hat

    hold_end <- min(reb_t + rebal_period - 1, T)
    port_ret_window <- returns[reb_t:hold_end, , drop=FALSE] %*% w_hat

    daily_port_ret <- c(daily_port_ret, port_ret_window)
    if (l == 1) {
      cum_rebal_returns[l] <- sum(port_ret_window)
    } else {
      cum_rebal_returns[l] <- cum_rebal_returns[l - 1] + sum(port_ret_window)
    }
  }


  N <- length(daily_port_ret)

  CER <- sum(daily_port_ret)

  mean_val <- CER / N
  devs <- daily_port_ret - mean_val
  stdev <- sqrt( sum(devs^2) / (N - 1) )

  SR <- (1/N) * (CER / stdev)

  list(
    rebal_dates              = rebalance_dates,
    weights                  = weights,
    daily_portfolio_returns  = daily_port_ret,
    cum_rebal_returns        = cum_rebal_returns,
    cumulative_excess_return = CER,
    standard_deviation       = stdev,
    sharpe_ratio             = SR
  )
}

#' @export
predict_portfolio <- function(
    returns,
    horizon = 1,
    bandwidth = 0.2,
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    lambda = 0.1
) {
  T <- nrow(returns)
  p <- ncol(returns)

  local_res <- localPCA(returns, bandwidth, max_factors, kernel_func)

  factor_cov <- cov(local_res$factors)

  residuals <- compute_residuals(local_res$factors, local_res$loadings, returns)
  residual_cov <- estimate_residual_cov(residuals, lambda)

  last_t <- T
  loadings_last <- local_res$loadings[[last_t]]
  Sigma_hat <- loadings_last %*% factor_cov %*% t(loadings_last) + residual_cov

  inv_cov <- solve(Sigma_hat)
  ones <- rep(1, p)
  w_unnorm <- inv_cov %*% ones
  w_hat <- as.numeric(w_unnorm / sum(w_unnorm))

  mean_returns <- colMeans(returns, na.rm = TRUE)
  expected_return <- sum(w_hat * mean_returns) * horizon

  risk <- sqrt(as.numeric(t(w_hat) %*% Sigma_hat %*% w_hat)) * sqrt(horizon)

  list(
    weights = w_hat,
    expected_return = expected_return,
    risk = risk
  )
}
