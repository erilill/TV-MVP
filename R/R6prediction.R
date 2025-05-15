# prediction

TVMVP$set("public", "predict_portfolio", function(
    horizon = 1, max_factors = 3,
    kernel_func = epanechnikov_kernel,
    min_return = NULL,
    max_SR = NULL,  # flag: if TRUE, compute maximum Sharpe portfolio
    rf = NULL
    ) {
  flag = TRUE
  if(is.null(private$data)) {
    cli::cli_alert_warning("data is empty")
    flag = FALSE
  }
  if(!flag) return(NULL) # return
  iT = private$iT; ip = private$ip

  # by default use the optimal m stored in the object, otherwise
  # determine optimal number of factors using Silverman’s bandwidth
  if(is.null(private$optimal_m)){
    m <- self$determine_factors(max_m=max_factors)
  } else {
    m <- private$optimal_m
  }
  # there should be bandwidth
  bandwidth <- private$bandwidth

  # Local PCA
  local_res <- localPCA(private$data, bandwidth, m, kernel_func)

  # Compute covariance matrix from local PCA results using POET
  Sigma_hat <- estimate_residual_cov_poet_local(localPCA_results = local_res,
                                                returns = private$data,
                                                M0 = 10,
                                                rho_grid = seq(0.005, 2, length.out = 30),
                                                floor_value = 1e-12,
                                                epsilon2 = 1e-6)$total_cov


  # Expected returns
  if (is.null(rf)) {
    mean_returns <- comp_expected_returns(private$data, horizon)
  } else {
    mean_returns <- comp_expected_returns(private$data, horizon) - rf
  }

  ## Compute MVP
  inv_cov <- chol2inv(chol(Sigma_hat))
  ones <- rep(1, ip)
  w_MVP_unnorm <- inv_cov %*% ones
  w_MVP <- as.numeric(w_MVP_unnorm / sum(w_MVP_unnorm))

  expected_return_MVP <- sum(w_MVP * mean_returns) * horizon
  risk_MVP <- sqrt(as.numeric(t(w_MVP) %*% Sigma_hat %*% w_MVP)) * sqrt(horizon)

  MVP <- list(
    weights = w_MVP,
    expected_return = expected_return_MVP,
    risk = risk_MVP,
    sharpe = (expected_return_MVP/horizon) / (risk_MVP/sqrt(horizon))
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
  method_names <- c("MVP")
  expected_returns_vec <- c(expected_return_MVP)
  risk_vec <- c(risk_MVP)
  sharpe_vec <- c(MVP$sharpe)

  if (!is.null(max_sr_portfolio)) {
    method_names <- c(method_names, "max_SR")
    expected_returns_vec <- c(expected_returns_vec, expected_return_max_sr)
    risk_vec <- c(risk_vec, risk_max_sr)
    sharpe_vec <- c(sharpe_vec, max_sr_portfolio$sharpe)
  }
  if (!is.null(constrained_portfolio)) {
    method_names <- c(method_names, "MVPConstrained")
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
    MVP = MVP,
    max_SR = if (!is.null(max_sr_portfolio)) max_sr_portfolio else NULL,
    MVPConstrained = if (!is.null(constrained_portfolio)) constrained_portfolio else NULL
  )

  return(out)
})
