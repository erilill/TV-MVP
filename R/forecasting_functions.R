#' @export
rolling_time_varying_mvp <- function(
    returns         ,
    initial_window  ,  # how many periods in the initial “estimation”
    rebal_period    ,  # holding window length (HT in the paper)
    max_factors     ,
    kernel_func    = epanechnikov_kernel,
    bandwidth_func = cv_bandwidth) {
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
      bandwidth <- handle_cv_bandwidth(returns, m, seq(0.05, 0.95, 0.05), kernel_func)
    }
    

    local_res <- localPCA(est_data, bandwidth, m, kernel_func)
    factor_cov   <- t(local_res$factors)%*%local_res$factors*(1/nrow(est_data))
    residuals    <- residuals(local_res$factors, local_res$loadings, est_data)
    residual_cov <- estimate_residual_cov(residuals)

    last_t <- nrow(est_data)
    loadings_mid <- local_res$loadings[[last_t]]

    Sigma_hat <- loadings_mid %*% factor_cov %*% t(loadings_mid) + residual_cov

    w_hat <- solve_minvar_portfolio(Sigma_hat)

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
#' @import quadprog
#' @export
predict_portfolio <- function(
    returns,
    horizon = 1,
    bandwidth_func = cv_bandwidth,
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    lambda = 0.1,
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
    bandwidth <- cv_bandwidth(est_data, m, seq(0.05, 0.95, 0.05), kernel_func)$optimal_h
    }

  # Perform Local PCA
  local_res <- localPCA(returns, bandwidth, m, kernel_func)

  # Estimate covariance matrices
  factor_cov <- t(local_res$factors) %*% local_res$factors / T
  res <- residuals(local_res$factors, local_res$loadings, returns)
  residual_cov <- estimate_residual_cov(res, lambda)

  # Compute Sigma_hat
  loadings_last <- local_res$loadings[[T]]
  Sigma_hat <- loadings_last %*% factor_cov %*% t(loadings_last) + residual_cov

  # Ensure Sigma_hat is well-conditioned
  Sigma_hat <- Sigma_hat + diag(1e-6, p)

  # **Corrected Mean Returns Calculation**
  mean_returns <- as.vector(t(loadings_last %*% local_res$factors[T,]))  # Ensures 1 x p

  ### **Global Minimum Variance Portfolio (GMVP)**
  inv_cov <- solve(Sigma_hat)
  ones <- rep(1, p)
  w_gmv_unnorm <- inv_cov %*% ones
  w_gmv <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))  # Normalize weights

  # Compute GMVP Expected Return and Risk
  expected_return_gmv <- sum(w_gmv * mean_returns) * horizon
  risk_gmv <- sqrt(as.numeric(t(w_gmv) %*% Sigma_hat %*% w_gmv)) * sqrt(horizon)

  ### **Minimum Variance Portfolio with Return Constraint**
  if (!is.null(min_return)) {
    Dmat <- as.matrix(Sigma_hat) + diag(1e-6, p)  # Regularization
    dvec <- rep(0, p)  # Corrected vector format
    Amat <- cbind(rep(1, p), mean_returns)  # Constraints matrix (p x 2)
    bvec <- c(1, min_return / horizon)  # Constraint values
    meq <- 2  # **Two equality constraints**

    result <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
    w_constrained <- as.numeric(result$solution)  # Extract optimized weights

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
#' @import quadprog
#' @export
solve_minvar_portfolio <- function(Sigma) {
  p <- ncol(Sigma)
  ones <- matrix(1, nrow = p, ncol = 1)

  # Solve using quadprog
  result <- solve.QP(
    Dmat = as.matrix(Sigma),
    dvec = rep(0, p),
    Amat = ones,
    bvec = 1,
    meq = 1
  )

  return(result$solution)  # Extract optimal weights
}

#' Helper function
find_smallest_matrix <- function(matrix_list) {
  if (length(matrix_list) == 0) {
  stop("The list is empty.")
    }
  # Extract dimensions of all matrices
  dims <- sapply(matrix_list, function(mat) c(nrow(mat), ncol(mat)))
  # Compute total number of elements (rows * cols)
  total_elements <- apply(dims, 2, prod)
  # Find the index of the smallest matrix
  smallest_index <- which.min(total_elements)
  # Return the smallest matrix
  return(dim(matrix_list[[smallest_index]]))
}

#' Helper function
handle_cv_bandwidth <- function(returns, m, candidate_h, kernel_func) {
  h <- try(cv_bandwidth(returns, m, candidate_h, kernel_func), silent = TRUE)
  
  if (inherits(h, "try-error")) {  # Check if an error occurred
    message("Error detected in bandwidth cross-validation. Attempting to lower m...")
    
    # Extract numeric values from the error message
    error_message <- as.character(h)
    m_eff_candidates <- as.numeric(unlist(regmatches(error_message, gregexpr("[0-9]+", error_message))))
    
    # Filter valid m_eff values (remove zeros, negatives, and NAs)
    valid_m_eff <- m_eff_candidates[m_eff_candidates > 0 & !is.na(m_eff_candidates)]
    
    if (length(valid_m_eff) > 0) {
      new_m <- min(valid_m_eff)  # Pick the smallest valid m_eff
      message(sprintf("Switching to m = %d.", new_m))
      
      return(cv_bandwidth(returns, new_m, candidate_h, kernel_func)$optimal_h)
    } else {
      stop("Failed to extract a valid m_eff from the error message. Bandwidth selection aborted.")
    }
  } 
  
  return(h$optimal_h)
}


