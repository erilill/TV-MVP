# expanding_tvmvp

TVMVP$set("public", "expanding_tvmvp", function(
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
  flag = TRUE
  if(is.null(private$data)) {
    cli::cli_alert_warning("data is empty")
    flag = FALSE
  }
  if(!flag) return(NULL) # return
  iT = private$iT; ip = private$ip

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
  # Determine when to update m_local based on 252-day intervals
  m_update_flags <- rep(FALSE, RT)
  days_since_last <- 0
  for (i in seq_len(RT)) {
    if (i == 1 || days_since_last >= 252) {
      m_update_flags[i] <- TRUE
      days_since_last <- 0
    } else {
      days_since_last <- days_since_last + rebal_period
    }
  }

  for (l in seq_len(RT)) {
    if (m_update_flags[l]){

      reb_t <- rebalance_dates[l]
      est_data <- private$data[1:(reb_t - 1), , drop = FALSE]
      # Use the estimation window to compute the bandwidth
      bandwidth <- silverman(est_data)
      m <- determine_factors(est_data, max_factors, bandwidth)$optimal_m
    }

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
    hold_len <- min(rebal_period, (iT - initial_window) - length(daily_port_ret_tvmvp))
    hold_end <- min(reb_t + hold_len - 1, iT)

    # Portfolio returns over the holding period:
    # TVMVP Portfolio returns
    port_ret_window_tvmvp <- as.numeric(private$data[reb_t:hold_end, , drop = FALSE] %*% w_hat)
    daily_port_ret_tvmvp <- c(daily_port_ret_tvmvp, port_ret_window_tvmvp)

    # Equal Weights Portfolio returns (baseline)
    w_equal <- rep(1/ip, ip)
    port_ret_window_equal <- as.numeric(private$data[reb_t:hold_end, , drop = FALSE] %*% w_equal)
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
    CER     = c(CER_tvmvp, CER_equal),
    MER     = c(mean_ret_tvmvp, mean_ret_equal),
    SD      = c(sd_tvmvp, sd_equal),
    SR      = c(SR_tvmvp, SR_equal),
    MER_ann = c(mean_annualized_tvmvp, mean_annualized_equal),
    SD_ann  = c(sd_annualized_tvmvp, sd_annualized_equal)
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
  out <- ExpandingWindow$new(
    summary = summary_df,
    TVMVP   = TVMVP,
    Equal   = Equal
  )
  class(out) <- c("ExpandingWindow", class(out))
  return(out)
})
