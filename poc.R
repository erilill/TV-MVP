# Load necessary libraries
library(MASS)        # For matrix operations
library(matrixcalc)  # For matrix calculations
library(glmnet)      # For lasso penalization
library(FactoMineR)  # For PCA
library(quadprog)    # For quadratic programming
library(forecast)    # For ARIMA comparison
library(zoo)         # For rolling operations
library(ggplot2)     # For visualization


################################################################################

# Example with simulated data
set.seed(123)
T <- 200  # Number of time periods
p <- 50   # Number of assets
returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

# Ensure data is properly scaled (optional)
returns <- scale(returns)

################################################################################
# Kernel function (Epanechnikov)
epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}

# Boundary kernel
boundary_kernel <- function(t, r, T, h, kernel_func) {
  scaled_diff <- (t - r) / (T*h)  # Argument of the kernel
  k_val <- kernel_func(scaled_diff) / h  # Scale by h
  
  # Determine the region of r
  Th_floor <- floor(T * h)
  
  if (r < Th_floor) {
    # Lower boundary case
    integral_val <- integrate(kernel_func, lower = -r / (T * h), upper = 1)$value # articles are conflicting, upper = 1 vs. Inf
    return(k_val / integral_val)
  } else if (r > (T - Th_floor)) {
    # Upper boundary case
    integral_val <- integrate(kernel_func, lower = -1, upper = (1 - r / T) / h)$value # articles are conflicting, lower = -1 vs. -Inf
    return(k_val / integral_val)
  } else {
    # Middle region
    return(k_val)
  }
}

# Two-fold convolution kernel
two_fold_convolution_kernel <- function(u, kernel_func) {
  result <- ifelse(
    abs(u) <= 2, # Not entirely sure about this, maybe it is dependent on kernel_func?
    sapply(u, function(u_val) {
      integrand <- function(v) kernel_func(v) * kernel_func(u_val - v)
      integrate(integrand, lower = -1, upper = 1)$value
    }),
    0
  )
  return(result)
}

# Function to Compute V_m (Sum of Squared Residuals)
compute_V_m <- function(returns, m, kernel_func, bandwidth) {
  p <- ncol(returns)
  T <- nrow(returns)
  total_ssr <- 0  # Sum of squared residuals
  
  for (x in 1:T) {
    # Compute boundary kernel weights for time x
    w_x <- sapply(1:T, function(t) boundary_kernel(x, t, T, bandwidth, kernel_func))
    w_x <- w_x / sum(w_x)  # Normalize weights
    
    # Apply weights to returns
    sqrt_w_x <- sqrt(w_x)
    weighted_returns <- sweep(returns, 1, sqrt_w_x, `*`)
    
    # Perform PCA
    pca_result <- prcomp(weighted_returns, center = FALSE, scale. = FALSE)
    num_pcs <- min(m, ncol(pca_result$x))
    if (num_pcs < 1) next  # Skip if no PCs are found
    
    # Extract factor scores and loadings
    Fhat <- pca_result$x[, 1:num_pcs, drop = FALSE] / sqrt(T)  # Normalize
    loadings_hat <- pca_result$rotation[, 1:num_pcs, drop = FALSE]
    
    # Compute fitted values and residuals
    fitted <- Fhat %*% t(loadings_hat)
    Resid_x <- returns - fitted
    total_ssr <- total_ssr + sum(Resid_x^2)
  }
  
  # Compute V_m
  V_m <- total_ssr / (p * T)
  return(V_m)
}

# Function to Select Optimal Number of Factors using Information Criterion
select_optimal_factors <- function(returns, max_factors, T_h, kernel_func, bandwidth) {
  p <- ncol(returns)
  T <- nrow(returns)
  
  IC_values <- numeric(max_factors)
  V_m_values <- numeric(max_factors)
  penalty_values <- numeric(max_factors)
  
  for (m in 1:max_factors) {
    V_m <- compute_V_m(returns, m, kernel_func, bandwidth)
    V_m_values[m] <- V_m
    
    # Compute penalty
    penalty <- (p + T_h) / (p * T_h) * log((p * T_h) / (p + T_h)) * m
    penalty_values[m] <- penalty
    
    # Compute Information Criterion (IC)
    IC_values[m] <- log(V_m) + penalty
  }
  
  optimal_m <- which.min(IC_values)
  
  return(list(optimal_m = optimal_m, IC_values = IC_values, V_m_values = V_m_values, penalty_values = penalty_values))
}


# Function to Perform Local PCA
local_pca <- function(returns, r, bandwidth, m, kernel_func) {
  T <- nrow(returns)
  
  # Compute boundary kernel weights for time r
  w_r <- sapply(1:T, function(t) boundary_kernel(r, t, T, bandwidth, kernel_func))
  w_r <- w_r / sum(w_r)  # Normalize weights
  
  # Apply weights to returns
  sqrt_w_r <- sqrt(w_r)
  weighted_returns <- sweep(returns, 1, sqrt_w_r, `*`)
  
  # Perform PCA
  pca_result <- prcomp(weighted_returns, center = FALSE, scale. = FALSE)
  
  # Determine actual number of factors
  num_factors <- min(m, ncol(pca_result$x))
  if (num_factors < 1) return(NULL)  # Return NULL if no factors are found
  
  # Extract factor scores and loadings
  Fhat_all <- pca_result$x[, 1:num_factors, drop = FALSE] / sqrt(T)  # Normalize
  loadings_all <- pca_result$rotation[, 1:num_factors, drop = FALSE]
  
  return(list(factors_full = Fhat_all, loadings_full = loadings_all))
}


# Apply local PCA over all time points
bandwidth <- (2.35/sqrt(12))*T^(-0.2)*p^(-0.1) # Silverman's rule of thumb
m <- select_optimal_factors(returns, max_factors = 10, T_h = T*bandwidth, kernel_func = epanechnikov_kernel, bandwidth)$optimal_m # Number of factors
factors_list <- list()
loadings_list <- list()

for (t in 1:T) {
  pca_result <- local_pca(returns, t, bandwidth, m, epanechnikov_kernel)
  factors_list[[t]] <- pca_result$factors
  loadings_list[[t]] <- pca_result$loadings
}

################################################################################
# Hypothesis testing

# A little unsure about the indexing of the local_loadings and local_factors here

compute_M_hat <- function(local_factors, global_factors, local_loadings, global_loadings, T, N, m) {
  M_hat <- 0
  if (m ==1){
    global_loadings <- matrix(global_loadings)
    global_factors <- matrix(global_factors)
  }
  for (i in 1:N) {
    for (t in 1:T) {
      common_H1 <- t(local_loadings[[t]][i, ]) %*% local_factors[[t]][t,]
      common_H0 <- (global_loadings[i, ]) %*% global_factors[t, ]
      M_hat <- M_hat + (common_H1 - common_H0)^2
    }
  }
  M_hat <- M_hat / (N * T)
  return(M_hat)
}


compute_B_pT <- function(local_factors, global_factors, residuals, h, T, p, kernel_func) {
  # 1) Precompute sum of residual squares per row s
  res2 <- rowSums(residuals^2)  # length T
  
  # 2) Kernel matrix K[s,t]
  K <- matrix(0, nrow=T, ncol=T)
  for (s in 1:T) {
    for (t in 1:T) {
      K[s, t] <- boundary_kernel(s, t, T, h, kernel_func)
    }
  }
  
  # 3) Local dot-product matrix, L[s,t] = l_s . l_t
  #    where l_s = local_factors[[s]][s, ], a length-r vector.
  L <- matrix(0, nrow=T, ncol=T)
  for (s in 1:T) {
    ls <- local_factors[[s]][s, ]
    for (t in 1:T) {
      lt <- local_factors[[t]][t, ]
      L[s, t] <- sum(ls * lt)
    }
  }
  
  # 4) Global dot-product matrix, G[s,t] = g_s . g_t
  #    if global_factors is T x r
  G <- global_factors %*% t(global_factors)
  
  # 5) Vectorized calculation of (K*L - G)^2, then multiply row s by res2[s], sum over all s,t
  D <- (K * L) - G  # elementwise
  D2 <- D^2
  val <- sum(D2 * res2[row(D2)])
  
  # 6) Final scaling
  B_pT <- (h^(1/2) / (T^2 * sqrt(p))) * val
  return(B_pT)
}


compute_V_pT <- function(local_factors, residuals, h, T, p, factor_cov, kernel_func) {
  V_pT <- 0
  for (s in 1:(T - 1)) {
    for (r in (s + 1):T) {
      k_bar_sr <- two_fold_convolution_kernel((s - r) / (T * h), kernel_func)
      term <- k_bar_sr^2 * (t(local_factors[[s]][s,]) %*% factor_cov %*% local_factors[[r]][r,])^2
      V_pT <- V_pT + term * (t(residuals[s, ]) %*% residuals[r, ])^2
    }
  }
  V_pT <- (2 / (T^2 * p * h)) * V_pT
  return(V_pT)
}


compute_J_pT <- function(B_pT, V_pT, M_hat, T, p, h) {
  J_pT <- (T * sqrt(p) * sqrt(h) * M_hat - B_pT) / sqrt(V_pT)
  return(J_pT)
}

################################################################################


# Factor covariance matrix
factor_covariance <- cov(do.call(rbind, factors_list))  # Aggregate over all factors


# Estimate residuals
residuals <- matrix(NA, nrow = T, ncol = p)
for (t in 1:T) {
  modeled_returns_t <- factors_list[[t]] %*% t(loadings_list[[t]])
  residuals[t, ] <- returns[t, ] - modeled_returns_t[t,]
}


# Function to Estimate Residual Covariance with Lasso Penalization
estimate_residual_cov <- function(residuals, lambda) {
  p <- ncol(residuals)
  
  # Compute the sample covariance matrix
  S <- cov(residuals)
  
  # Initialize the sparse covariance matrix
  sparse_cov <- matrix(0, nrow = p, ncol = p)
  
  # Apply Lasso regression row-wise
  for (i in 1:p) {
    response <- S[i, -i]
    predictors <- S[-i, -i]
    
    # Fit Lasso regression without intercept
    fit <- glmnet(predictors, response, alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE)
    
    # Extract coefficients and assign to sparse covariance matrix
    coef_i <- as.vector(coef(fit, s = lambda))[-1]  # Exclude intercept
    sparse_cov[i, -i] <- coef_i
  }
  
  # Symmetrize the covariance matrix
  sparse_cov <- (sparse_cov + t(sparse_cov)) / 2
  
  # Assign original diagonal elements
  diag(sparse_cov) <- diag(S)
  
  return(sparse_cov)
}

lambda <- 0.1  # Set penalty parameter
residual_covariance <- estimate_residual_cov(residuals, lambda)


# Compute Time-Varying Covariance Matrices
compute_time_varying_cov <- function(loadings, factor_cov, residual_cov) {
  return(loadings %*% factor_cov %*% t(loadings) + residual_cov)
}

time_varying_cov_list <- mapply(
  compute_time_varying_cov,
  loadings = loadings_list,
  MoreArgs = list(factor_cov = factor_covariance, residual_cov = residual_covariance),
  SIMPLIFY = FALSE
)



compute_optimal_weights <- function(cov_matrix, p) {
  if (any(is.na(cov_matrix)) || !is.positive.definite(round(cov_matrix, 10))) {
    # Add a small ridge term to stabilize inversion
    cov_matrix <- cov_matrix + diag(1e-6, p)
  }
  
  Dmat <- solve(cov_matrix)
  dvec <- rep(0, p)
  Amat <- matrix(1, nrow = p, ncol = 1)  # Constraint: sum(weights) = 1
  bvec <- 1
  
  # Solve QP
  result <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
    error = function(e) NULL
  )
  
  if (!is.null(result)) {
    return(result$solution)
  } else {
    # Return equal weights if QP fails
    return(rep(1/p, p))
  }
}
optimal_weights <- lapply(time_varying_cov_list, function(cov_mat) {
  compute_optimal_weights(cov_mat, p)
})


compute_portfolio_risk <- function(weights, cov_matrix) {
  return(sqrt(as.numeric(t(weights) %*% cov_matrix %*% weights)))
}
portfolio_risk <- mapply(
  compute_portfolio_risk,
  weights = optimal_weights,
  cov_matrix = time_varying_cov_list,
  SIMPLIFY = TRUE
)



# Compute Expected Returns (Mean of Returns)
expected_returns <- colMeans(returns, na.rm = TRUE)

# Compute Portfolio Sharpe Ratios
compute_sharpe_ratio <- function(weights, expected_returns, cov_matrix) {
  portfolio_return <- sum(weights * expected_returns, na.rm = TRUE)
  portfolio_std <- sqrt(as.numeric(t(weights) %*% cov_matrix %*% weights))
  if (portfolio_std == 0) return(NA)
  return(portfolio_return / portfolio_std)
}

portfolio_sharpe <- mapply(
  compute_sharpe_ratio,
  weights = optimal_weights,
  cov_matrix = time_varying_cov_list,
  MoreArgs = list(expected_returns = expected_returns),
  SIMPLIFY = TRUE
)




################################################################################
# Testing hypothesis test

# Global PCA
global_pca <- prcomp(returns, scale. = FALSE, center = FALSE)
global_factors <- global_pca$x[, 1:m]/sqrt(T)
global_loadings <- global_pca$rotation[, 1:m]

# Compute M_hat
M_hat <- compute_M_hat(factors_list, global_factors, loadings_list, global_loadings, T, p, m)
print(paste("M_hat:", M_hat)) #"M_hat: 0.00019855543777811"

# Compute B_pT
B_pT <- compute_B_pT(factors_list, global_factors, residuals, bandwidth, T, p, epanechnikov_kernel)
print(paste("B_pT:", B_pT)) # 0.000326884513912354

# Compute V_pT
V_pT <- compute_V_pT(factors_list, residuals, bandwidth, T, p, factor_covariance, epanechnikov_kernel)
print(paste("V_pT:", V_pT)) # V_pT: 2.82995713688749e-15

# Compute J_pT
J_pT <- compute_J_pT(B_pT, V_pT, M_hat, T, p, bandwidth)
print(paste("J_pT:", J_pT)) # J_pT: 2098573.53331389

################################################################################
# Forecasting
forecast_local_factor_model <- function(
    returns,
    W,
    m,
    bandwidth,
    kernel_func = epanechnikov_kernel,
    lambda = 0.1  # Penalization parameter
) {
  T <- nrow(returns)
  p <- ncol(returns)
  
  # Initialize storage
  forecasts <- matrix(NA, nrow = T, ncol = p)
  est_covariances <- vector("list", T)
  residual_covariances <- vector("list", T)
  optimal_weights <- vector("list", T)
  
  for (tau in (W + 1):T) {
    # Compute boundary kernel weights for times 1..(tau-1)
    weight_vec <- sapply(1:(tau - 1), function(r) {
      boundary_kernel(tau, r, tau - 1, bandwidth, kernel_func)
    })
    weight_vec <- weight_vec / sum(weight_vec)  # Normalize weights
    sqrt_w <- sqrt(weight_vec)
    
    # Apply weights to sub-returns
    weighted_subreturns <- sweep(returns[1:(tau - 1), ], 1, sqrt_w, `*`)
    
    # Perform Local PCA
    pca_result <- prcomp(weighted_subreturns, center = FALSE, scale. = FALSE)
    num_factors <- min(m, ncol(pca_result$x))
    if (num_factors < 1) next  # Skip if no factors
    
    # Extract factors and loadings
    Fhat <- pca_result$x[, 1:num_factors, drop = FALSE] / sqrt(tau - 1)  # Normalize
    loadings_hat <- pca_result$rotation[, 1:num_factors, drop = FALSE]
    
    # Forecast for time tau: Use last factor scores
    F_tau_hat <- Fhat[nrow(Fhat), , drop = FALSE]
    R_hat_tau <- loadings_hat %*% t(F_tau_hat)
    forecasts[tau, ] <- as.vector(R_hat_tau)
    
    # Estimate Factor Covariance
    factor_cov <- cov(Fhat, use = "complete.obs")
    
    # Estimate Residual Covariance
    residuals_sub <- returns[1:(tau - 1), ] - Fhat %*% t(loadings_hat)
    resid_cov <- estimate_residual_cov(residuals_sub, lambda)
    residual_covariances[[tau]] <- resid_cov
    
    # Combine to get Estimated Covariance
    est_cov <- loadings_hat %*% factor_cov %*% t(loadings_hat) + resid_cov
    est_covariances[[tau]] <- est_cov
    
    # Compute Optimal Weights using QP
    weights <- compute_optimal_weights(est_cov, p)
    optimal_weights[[tau]] <- weights
  }
  
  return(list(
    forecasts = forecasts,
    est_covariances = est_covariances,
    residual_covariances = residual_covariances,
    optimal_weights = optimal_weights
  ))
}



# Evaluation Function
evaluate_forecasts <- function(
    returns, 
    forecasts,                  # T x p matrix of predicted means
    est_covariances,            # list of length T, each a p x p predicted covariance
    residual_covariances_pred,  # list of length T, each a p x p predicted residual covariance
    weights_est,                # list of length T, each a p-vector of portfolio weights
    window_eval,                # Vector of time indices, e.g., (W +1):T
    realized_covariances   = NULL,       # list of length T
    realized_resid_covariances = NULL,   # list of length T
    true_weights           = NULL,       # list of length T
    realized_sharpes       = NULL        # numeric vector of length T
) {
  # Initialize metrics
  metrics <- list(
    risk_error = NA,
    weight_error = NA,
    sharpe_error = NA,
    covariance_error = NA,
    residual_cov_error = NA
  )
  
  # 1. Risk Error
  if (!is.null(realized_covariances) && !is.null(weights_est)) {
    risk_diffs <- mapply(function(w, est_cov, real_cov) {
      if (is.null(w) || is.null(est_cov) || is.null(real_cov)) return(NA)
      pred_risk <- sqrt(as.numeric(t(w) %*% est_cov %*% w))
      real_risk <- sqrt(as.numeric(t(w) %*% real_cov %*% w))
      return(abs(real_risk - pred_risk))
    }, 
    w = weights_est[window_eval], 
    est_cov = est_covariances[window_eval], 
    real_cov = realized_covariances[window_eval])
    
    metrics$risk_error <- mean(risk_diffs, na.rm = TRUE)
  }
  
  # 2. Weight Error
  if (!is.null(true_weights) && !is.null(weights_est)) {
    weight_diffs <- mapply(function(w_est, w_true) {
      if (is.null(w_est) || is.null(w_true)) return(NA)
      return(sqrt(sum((w_est - w_true)^2)))
    },
    w_est = weights_est[window_eval],
    w_true = true_weights[window_eval])
    
    metrics$weight_error <- mean(weight_diffs, na.rm = TRUE)
  }
  
  # 3. Sharpe Ratio Error
  if (!is.null(realized_sharpes) && !is.null(weights_est) && !is.null(forecasts)) {
    sharpe_diffs <- mapply(function(w, forecast_row, realized_sr, est_cov) {
      if (is.null(w) || is.null(forecast_row) || is.na(realized_sr) || is.null(est_cov)) return(NA)
      mu_hat <- sum(forecast_row * w, na.rm = TRUE)
      sigma_hat <- sqrt(as.numeric(t(w) %*% est_cov %*% w))
      if (sigma_hat == 0 || is.na(sigma_hat)) return(NA)
      sr_hat <- mu_hat / sigma_hat
      return(abs(realized_sr - sr_hat))
    },
    w = weights_est[window_eval],
    forecast_row = split(forecasts[window_eval, ], row(forecasts[window_eval, ])),
    realized_sr = realized_sharpes[window_eval],
    est_cov = est_covariances[window_eval])
    
    metrics$sharpe_error <- mean(sharpe_diffs, na.rm = TRUE)
  }
  
  # 4. Covariance Error
  if (!is.null(realized_covariances) && !is.null(est_covariances)) {
    cov_diffs <- mapply(function(est_cov, real_cov) {
      if (is.null(est_cov) || is.null(real_cov)) return(NA)
      return(norm(est_cov - real_cov, type = "F"))
    },
    est_cov = est_covariances[window_eval],
    real_cov = realized_covariances[window_eval])
    
    metrics$covariance_error <- mean(cov_diffs, na.rm = TRUE)
  }
  
  # 5. Residual Covariance Error
  if (!is.null(realized_resid_covariances) && !is.null(residual_covariances_pred)) {
    resid_cov_diffs <- mapply(function(pred_resid, real_resid) {
      if (is.null(pred_resid) || is.null(real_resid)) return(NA)
      return(norm(pred_resid - real_resid, type = "F"))
    },
    pred_resid = residual_covariances_pred[window_eval],
    real_resid = realized_resid_covariances[window_eval])
    
    metrics$residual_cov_error <- mean(resid_cov_diffs, na.rm = TRUE)
  }
  
  return(metrics)
}

W <- 100
bandwidth <- (2.35 / sqrt(12)) * T^(-0.2) * p^(-0.1)

# 3) Get forecasts and estimated covariances from your local factor model
oos_results <- forecast_local_factor_model(
  returns      = returns,
  W            = W,
  m            = m,
  bandwidth    = bandwidth,
  kernel_func  = epanechnikov_kernel
)

# 4) Realized Covariances
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

# 5) Optimal weights
optimal_weights <- list()
for (t in (W+1):T) {
  time_varying_cov <- loadings_list[[t]] %*% factor_covariance %*% t(loadings_list[[t]]) + residual_covariance
  Dmat <- solve(time_varying_cov)  # Inverse covariance matrix
  dvec <- rep(0, p)
  Amat <- cbind(rep(1, p))  # Constraint: sum of weights = 1
  bvec <- 1
  
  # Solve quadratic programming problem
  result <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  optimal_weights[[t]] <- result$solution
}

# 6) Sharpe ratio
realized_sharpes <- rep(NA, T)
for (t in (W+1):T) {
  tmp <- returns[(t - W):t, , drop = FALSE]
  
  # Weighted portfolio (assuming you already have weights and time-varying covariance)
  weights <- optimal_weights[[t]]
  time_varying_cov <- loadings_list[[t]] %*% factor_covariance %*% t(loadings_list[[t]]) + residual_covariance
  
  # Realized returns of the portfolio from day t-W..t
  rets_t <- rowSums( sweep(tmp, 2, weights, `*`) )  # a vector
  
  # Predicted or "current" portfolio std dev using the factor-based covariance
  rets_t_std <- sqrt( t(weights) %*% time_varying_cov %*% weights )
  
  realized_sharpes[t] <- mean(rets_t) / rets_t_std
}



# 7) Evaluate
perf <- evaluate_forecasts(
  returns                      = returns,
  forecasts                    = oos_results$forecasts,
  est_covariances              = oos_results$est_covariances,
  residual_covariances_pred    = oos_results$residual_covariances,
  weights_est                  = oos_results$optimal_weights, 
  window_eval                  = (W +1):T,
  realized_covariances         = realized_covariances,
  realized_resid_covariances   = realized_resid_covariances,
  true_weights                 = optimal_weights,
  realized_sharpes             = realized_sharpes
)

print(perf)

#$risk_error
#[1] 0.1072943

#$weight_error
#[1] 0.04058218

#$sharpe_error
#[1] 0.1262457

#$covariance_error
#[1] 12.84417

#$residual_cov_error
#[1] 12.78617
################################################################################
# Testing functions

# Test boundary kernel weights
test_boundary_kernel <- function() {
  T <- 200
  h <- 0.1
  weights <- sapply(1:T, function(t) boundary_kernel(t, T, T, h, epanechnikov_kernel))
  weights <- weights/sum(weights)
  
  # Check if weights sum approximately to 1
  sum_weights <- sum(weights)
  print(paste("Sum of weights (should be close to 1):", sum_weights))
  
  # Plot weights for visual inspection
  plot(1:T, weights, type = "l", main = "Boundary Kernel Weights", xlab = "Time", ylab = "Weight")
}

# Test local PCA function
test_local_pca <- function() {
  set.seed(123)
  T <- 200
  p <- 50
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
  bandwidth <- 0.1
  
  pca_result <- local_pca(returns, r = T, bandwidth = bandwidth, m = m, kernel_func = epanechnikov_kernel)
  
  # Print dimensions of factors and loadings
  print(paste("Factors dimension (should be m):", length(pca_result$factors)))
  print(paste("Loadings dimension (should be p x m):", dim(pca_result$loadings)[1], "x", dim(pca_result$loadings)[2]))
}


# Test information criterion function
test_information_criterion <- function() {
  set.seed(123)
  T <- 200
  p <- 50
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
  
  result <- select_optimal_factors(returns, max_factors = 10, T_h = T*bandwidth, epanechnikov_kernel ,bandwidth)
  
  # Print results
  print(paste("Optimal number of factors:", result$optimal_m))
  print("IC values:")
  print(result$IC_values)
  
  # Plot IC values
  plot(1:10, result$IC_values, type = "b", main = "Information Criterion vs. Number of Factors",
       xlab = "Number of Factors", ylab = "IC")
}




# Test hypothesis testing components
test_hypothesis_testing <- function() {
  set.seed(123)
  T <- 200
  p <- 50
  m <- 3  # Number of factors (adjust as needed)
  bandwidth <- (2.35 / sqrt(12)) * T^(-0.2) * p^(-0.1)  # Silverman's rule of thumb
  
  # Simulate returns
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
  
  # Generate mock local factors and loadings (adjusted to new setup)
  factors_list <- lapply(1:T, function(t) matrix(rnorm(T * m), nrow = T, ncol = m))
  loadings_list <- lapply(1:T, function(t) matrix(rnorm(p * m), nrow = p, ncol = m))
  
  # Generate residuals
  residuals <- matrix(rnorm(T * p), nrow = T, ncol = p)
  
  # Simulate global factors and loadings
  global_factors <- matrix(rnorm(T * m), nrow = T, ncol = m)
  global_loadings <- matrix(rnorm(p * m), nrow = p, ncol = m)
  
  # Mock factor covariance matrix
  factor_covariance <- diag(m)
  
  # Compute hypothesis testing components
  M_hat <- compute_M_hat(factors_list, global_factors, loadings_list, global_loadings, T, p, m)
  B_pT <- compute_B_pT(factors_list, global_factors, residuals, bandwidth, T, p, epanechnikov_kernel)
  V_pT <- compute_V_pT(factors_list, residuals, bandwidth, T, p, factor_covariance, epanechnikov_kernel)
  J_pT <- compute_J_pT(B_pT, V_pT, M_hat, T, p, bandwidth)
  
  # Print results
  print(paste("M_hat:", M_hat))
  print(paste("B_pT:", B_pT))
  print(paste("V_pT:", V_pT))
  print(paste("J_pT:", J_pT))
}



test_forecasting <- function() {
  set.seed(123)
  T <- 200
  p <- 50
  W <- 100  # Window length
  m <- 3    # Number of factors
  bandwidth <- (2.35 / sqrt(12)) * T^(-0.2) * p^(-0.1)  # Silverman's rule of thumb
  
  # Simulated returns data
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
  
  # Initialize matrix to store forecasts
  forecasts <- matrix(NA, nrow = T, ncol = p)
  
  for (tau in (W + 1):T) {
    # Calculate weights using the boundary kernel
    weight_vec <- sapply(1:(tau - 1), function(r) boundary_kernel(tau - 1, r, tau - 1, bandwidth, epanechnikov_kernel))
    weight_vec <- weight_vec / sum(weight_vec)  # Normalize weights
    sqrt_w <- sqrt(weight_vec)
    
    # Weighted returns up to tau - 1
    weighted_subreturns <- sweep(returns[1:(tau - 1), , drop = FALSE], 1, sqrt_w, `*`)
    
    # Perform local PCA
    pca_result <- prcomp(weighted_subreturns, center = FALSE, scale = FALSE)
    
    # Check the number of components returned by PCA
    num_factors <- min(m, ncol(pca_result$x))  # Ensure we don't exceed available factors
    if (num_factors == 0) {
      warning(paste("Skipping tau =", tau, "due to insufficient factors."))
      next
    }
    
    # Extract factors and loadings
    local_factors_mat <- pca_result$x[, 1:num_factors, drop = FALSE]  # (tau-1) x num_factors
    loadings_hat <- t(local_factors_mat) %*% weighted_subreturns / (tau - 1)  # num_factors x p
    loadings_hat <- t(loadings_hat)  # p x num_factors
    
    # Get the last factor row (F_tau_hat)
    if (num_factors == 1) {
      F_tau_hat <- local_factors_mat[nrow(local_factors_mat), drop = FALSE]  # Ensure it's a matrix
    } else {
      F_tau_hat <- local_factors_mat[nrow(local_factors_mat), ]
    }
    
    # Forecast returns at time tau
    Rhat_tau <- loadings_hat %*% F_tau_hat
    
    # Store forecast for time tau
    forecasts[tau, ] <- Rhat_tau
  }
  
  # Sanity check: Forecast dimensions and sample forecast
  print(dim(forecasts))  # Should be T x p
  print(forecasts[(W + 1):(W + 5), 1:5])  # Print a few forecasts for the first 5 assets
  
  # Optionally, return the forecasts for further analysis
  return(forecasts)
}


validate_convolution_kernel <- function(kernel_func) {
  integrand <- function(u) two_fold_convolution_kernel(u, kernel_func)
  result <- integrate(integrand, lower = -2, upper = 2)
  print(paste("Integral of two-fold convolution kernel:", result$value))
}

validate_convolution_kernel(epanechnikov_kernel)
test_boundary_kernel()
test_local_pca()
test_information_criterion()
test_hypothesis_testing()
test_forecasting()

