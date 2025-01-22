library(MASS)       # For matrix operations
library(matrixcalc) # For matrix calculations
library(glmnet)     # For lasso penalization
library(FactoMineR) # For PCA
library(quadprog)   # For quadratic programming
library(FE)         # Testing with real data

################################################################################

# Example with simulated data
set.seed(123)
T <- 200  # Number of time periods
p <- 50   # Number of assets
returns <- portfolio_m[,5:124] #matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
T <- nrow(returns)
p <- ncol(returns)

# Ensure data is properly scaled (optional)
returns <- scale(returns)

################################################################################

# Kernel function (Epanechnikov)
epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}

boundary_kernel <- function(t, r, T, h, kernel_func) {
  scaled_diff <- (t - r) / (T*h)  # Argument of the kernel
  k_val <- kernel_func(scaled_diff) / h  # Scale by h
  
  # Determine the region of r
  Th_floor <- floor(T * h)
  
  if (r < Th_floor) {
    # Lower boundary case
    integral_val <- integrate(kernel_func, lower = -r / (T * h), upper = 1)$value
    return(k_val / integral_val)
    } else if (r > (T - Th_floor)) {
    # Upper boundary case
    integral_val <- integrate(kernel_func, lower = -1, upper = (1 - r / T) / h)$value
    return(k_val / integral_val)
      } else {
    # Middle region
    return(k_val)
  }
}

# Two-fold convolution kernel
two_fold_convolution_kernel <- function(u, kernel_func) {
  result <- ifelse(
    abs(u) <= 2,
    sapply(u, function(u_val) {
      integrand <- function(v) kernel_func(v) * kernel_func(u_val - v)
      integrate(integrand, lower = -1, upper = 1)$value
    }),
    0
  )
  return(result)
}

# Function to estimate IC(m) and find the optimal number of factors
compute_V_m <- function(returns, m, kernel_func, bandwidth, ridge = 1e-6) {
  p <- ncol(returns)
  T <- nrow(returns)
  
  # Compute weights for each time point
  weights <- matrix(0, nrow = T, ncol = T)
  for (x in 1:T) {
    weights[, x] <- sapply(1:T, function(t) boundary_kernel(t, x, T, bandwidth, kernel_func))
  }
  
  # Initialize residual sum of squares
  residual_sum <- 0
  
  for (x in 1:T) {
    # Weight the returns
    sqrt_weights <- sqrt(weights[, x])
    weighted_returns <- sweep(returns, 1, sqrt_weights, `*`)
    
    # PCA on weighted returns
    pca_result <- prcomp(weighted_returns, scale. = FALSE, center = FALSE)
    factors <- pca_result$x[, 1:m]  # Dimensions: T x m if m > 1
    loadings <- pca_result$rotation[, 1:m]  # Dimensions: p x m
    
    # Handle case where m = 1
    if (m == 1) {
      factors <- matrix(factors, ncol = 1)  # Convert vector to matrix
      loadings <- matrix(loadings, ncol = 1)  # Convert vector to matrix
    }
    
    # Optimal loadings for time x
    R_x <- returns[x, , drop = FALSE]  # Dimensions: 1 x p
    F_x <- factors[x, , drop = FALSE]  # Dimensions: 1 x m
    
    # Regularize t(F_x) %*% F_x
    F_x_t_F_x <- t(F_x) %*% F_x + ridge * diag(m)
    
    # Compute B_x with regularization
    B_x <- solve(F_x_t_F_x) %*% t(F_x) %*% R_x  # Dimensions: m x 1
    
    # Compute residuals
    residuals <- R_x - F_x %*% B_x
    residual_sum <- residual_sum + sum(residuals^2)
  }
  
  # Final V_m
  V_m <- residual_sum / (p * T)
  return(V_m)
}


select_optimal_factors <- function(returns, max_factors, T_h, kernel_func, bandwidth) {
  p <- ncol(returns)
  T <- nrow(returns)
  
  IC_values <- numeric(max_factors)
  V_m_values <- numeric(max_factors)
  penalty_values <- numeric(max_factors)
  
  for (m in 1:max_factors) {
    # Compute V_m using the updated function
    V_m <- compute_V_m(returns, m, kernel_func, bandwidth)
    V_m_values[m] <- V_m
    
    # Compute penalty
    penalty <- (p + T_h) / (p * T_h) * log((p * T_h) / (p + T_h)) * m
    penalty_values[m] <- penalty
    
    # Compute IC
    IC_values[m] <- log(V_m) + penalty
  }
  
  optimal_m <- which.min(IC_values)
  
  #print(data.frame(m = 1:max_factors, V_m = V_m_values, Penalty = penalty_values, IC = IC_values, optimal_m = optimal_m))
  return(list(optimal_m = optimal_m, IC_values = IC_values))
}

# Local PCA function
local_pca <- function(returns, r, bandwidth, m, kernel_func) {
  T <- nrow(returns) 
  weights <- sapply(1:T, function(t) boundary_kernel(r, t, T, bandwidth, kernel_func))
  weights <- weights/sum(weights)
  sqrt_weights <- sqrt(weights)
  weighted_returns <- sweep(returns, 1, sqrt_weights, `*`)
  pca_result <- prcomp(weighted_returns, scale. = FALSE, center = FALSE) #Scale?
  factors <- pca_result$x[, 1:m]/sqrt(T)
  loadings <- t(factors) %*% weighted_returns / T
  
  # Handle case where m = 1
  if (m == 1) {
    factors <- matrix(factors, ncol = 1)  # Convert vector to matrix
    loadings <- matrix(loadings, ncol = 1)  # Convert vector to matrix
  }
  
  list(factors = factors[r,], loadings = t(loadings)) # correct? should I include all factors?
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

compute_M_hat <- function(local_factors, global_factors, local_loadings, global_loadings, T, N) {
  M_hat <- 0
  for (i in 1:N) {
    for (t in 1:T) {
      common_H1 <- local_loadings[[t]][i, ] %*% local_factors[t,]
      common_H0 <- global_loadings[i, ] %*% global_factors[t, ]
      M_hat <- M_hat + (common_H1 - common_H0)^2
    }
  }
  M_hat <- M_hat / (N * T)
  return(M_hat)
}


compute_B_pT <- function(local_factors, global_factors, residuals, h, T, p, kernel_func) {
  B_pT <- 0
  for (i in 1:p) {
    for (t in 1:T) {
      for (s in 1:T) {
        k_h_st <- boundary_kernel(s, t, T, h, kernel_func)
        diff <- (k_h_st * ((local_factors[s,]) %*% local_factors[t,]) - # Outer poroduct or dot product?
                   ((global_factors[s, ]) %*% global_factors[t, ]))^2
        B_pT <- B_pT + diff * residuals[s, i]^2
      }
    }
  }
  B_pT <- (h^(1/2) / (T^2 * sqrt(p))) * B_pT
  return(B_pT)
}

compute_V_pT <- function(local_factors, residuals, h, T, p, factor_cov, kernel_func) {
  V_pT <- 0
  for (s in 1:(T - 1)) {
    for (r in (s + 1):T) {
      k_bar_sr <- two_fold_convolution_kernel((s - r) / (T * h), kernel_func)
      term <- k_bar_sr^2 * (t(local_factors[s,]) %*% factor_cov %*% local_factors[r,])^2
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
  modeled_returns_t <- loadings_list[[t]] %*% matrix(factors_list[[t]], ncol = 1)
  residuals[t, ] <- returns[t, ] - as.vector(modeled_returns_t)
}


# Penalized covariance estimation for residuals
estimate_residual_cov <- function(residuals, lambda) {
  # Compute the sample covariance matrix
  S <- cov(residuals)
  p <- ncol(S)
  
  # Initialize the sparse covariance matrix
  sparse_cov <- matrix(0, nrow = p, ncol = p)
  
  # Apply Lasso penalization row by row
  # Is row by row correct?
  for (i in 1:p) {
    # Response variable: row i of the covariance matrix
    response <- S[i, -i]
    
    # Predictors: all other rows/columns except the diagonal
    predictors <- S[-i, -i]
    
    # Fit Lasso regression
    fit <- glmnet(predictors, response, alpha = 1, lambda = lambda)
    
    # Update sparse covariance matrix
    sparse_cov[i, -i] <- as.vector(coef(fit, s = lambda)[-1])  # Exclude intercept
  }
  
  # Symmetrize the covariance matrix
  sparse_cov <- (sparse_cov + t(sparse_cov)) / 2
  
  # Add diagonal elements from the original covariance
  diag(sparse_cov) <- diag(S)
  
  sparse_cov
}


lambda <- 0.1  # Set penalty parameter
residual_covariance <- estimate_residual_cov(residuals, lambda)


time_varying_cov <- list()
for (t in 1:T) {
  time_varying_cov[[t]] <- loadings_list[[t]] %*% factor_covariance %*% t(loadings_list[[t]]) + residual_covariance
}


optimal_weights <- list()
for (t in 1:T) {
  Dmat <- solve(time_varying_cov[[t]])  # Inverse covariance matrix
  dvec <- rep(0, p)
  Amat <- cbind(rep(1, p))  # Constraint: sum of weights = 1
  bvec <- 1
  
  # Solve quadratic programming problem
  result <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  optimal_weights[[t]] <- result$solution
}


portfolio_risk <- sapply(1:T, function(t) {
  weights <- optimal_weights[[t]]
  sqrt(t(weights) %*% time_varying_cov[[t]] %*% weights)
})


expected_returns <- apply(returns, 2, mean)  # Approximate expected returns
sharpe_ratios <- sapply(1:T, function(t) {
  weights <- optimal_weights[[t]]
  portfolio_return <- sum(weights * expected_returns)
  portfolio_std <- sqrt(t(weights) %*% time_varying_cov[[t]] %*% weights)
  portfolio_return / portfolio_std
})



################################################################################
# Testing hypothesis test

# Global PCA
global_pca <- prcomp(returns, scale. = FALSE, center = FALSE)
global_factors <- global_pca$x[, 1:m]/sqrt(T)
global_loadings <- global_pca$rotation[, 1:m]

# Compute M_hat
local_factors <- matrix(unlist(factors_list), ncol=m, byrow=T)
local_loadings_list <- lapply(1:T, function(t) loadings_list[[t]])

M_hat <- compute_M_hat(local_factors, global_factors, local_loadings_list, global_loadings, T, p)
print(paste("M_hat:", M_hat))

# Compute B_pT
B_pT <- compute_B_pT(local_factors, global_factors, residuals, bandwidth, T, p, epanechnikov_kernel)
print(paste("B_pT:", B_pT))

# Compute V_pT
V_pT <- compute_V_pT(local_factors, residuals, bandwidth, T, p, factor_covariance, epanechnikov_kernel)
print(paste("V_pT:", V_pT))

# Compute J_pT
J_pT <- compute_J_pT(B_pT, V_pT, M_hat, T, p, bandwidth)
print(paste("J_pT:", J_pT))

################################################################################
# Forecasting

# Suppose we want to forecast from t = T0+1 up to T.
# We'll do a loop from 'tau = T0+1' to 'T' and forecast 'returns[tau, ]'.

W <- 100  # window length for local smoothing, or might be entire [1..(tau-1)]

forecasts <- matrix(NA, nrow = T, ncol = p)

for (tau in (W+1):T) {
  weight_vec <- sapply(1:(tau-1), function(r) {
    boundary_kernel(tau-1, r, (tau-1), bandwidth, epanechnikov_kernel)
  })
  weight_vec <- weight_vec/sum(weight_vec)
  
  # Weighted returns: up to time tau-1
  sqrt_w <- sqrt(weight_vec)
  weighted_subreturns <- sweep(returns[1:(tau-1), , drop=FALSE], 1, sqrt_w, `*`)
  
  # Then do prcomp to get local factors:
  pca_result <- prcomp(weighted_subreturns, center=FALSE, scale=FALSE)
  
  local_factors_mat <- pca_result$x[, 1:m] # (tau-1) x m
  loadings_hat <- t(local_factors_mat) %*% weighted_subreturns / (tau-1)
  loadings_hat <- t(loadings_hat)  # p x m
  
  # For a naive approach:
  F_tau_hat <- local_factors_mat[nrow(local_factors_mat), ]  # last factor
  
  # forecast:
  Rhat_tau <- loadings_hat %*% F_tau_hat
  
  # store forecast for time tau
  forecasts[tau, ] <- Rhat_tau
}

# Initialize matrices for alternative forecasts
naive_mean_forecasts <- matrix(NA, nrow = T, ncol = p)
random_walk_forecasts <- matrix(NA, nrow = T, ncol = p)
ar1_forecasts <- matrix(NA, nrow = T, ncol = p)

# Loop over assets
for (j in 1:p) {
  # Historical Mean Forecast
  for (tau in (W+1):T) {
    naive_mean_forecasts[tau, j] <- mean(returns[1:(tau-1), j], na.rm = TRUE)
  }
  
  # Random Walk Forecast
  for (tau in (W+1):T) {
    random_walk_forecasts[tau, j] <- returns[tau-1, j]
  }
  
  # AR(1) Model Forecast
  for (tau in (W+1):T) {
    if ((tau - 1) > 1) {
      fit_ar1 <- arima(returns[1:(tau-1), j], order = c(1, 0, 0), include.mean = TRUE)
      ar1_forecasts[tau, j] <- predict(fit_ar1, n.ahead = 1)$pred
    }
  }
}

# Evaluate out-of-sample error for each method
actuals <- returns[(W+1):T, ]

# Your method
your_preds <- forecasts[(W+1):T, ]
your_errors <- actuals - your_preds
your_mse <- mean(your_errors^2, na.rm = TRUE)

# Historical Mean
naive_mean_preds <- naive_mean_forecasts[(W+1):T, ]
naive_mean_errors <- actuals - naive_mean_preds
naive_mean_mse <- mean(naive_mean_errors^2, na.rm = TRUE)

# Random Walk
random_walk_preds <- random_walk_forecasts[(W+1):T, ]
random_walk_errors <- actuals - random_walk_preds
random_walk_mse <- mean(random_walk_errors^2, na.rm = TRUE)

# AR(1)
ar1_preds <- ar1_forecasts[(W+1):T, ]
ar1_errors <- actuals - ar1_preds
ar1_mse <- mean(ar1_errors^2, na.rm = TRUE)

# Compare Results
mse_results <- data.frame(
  Method = c("Your Method", "Historical Mean", "Random Walk", "AR(1)"),
  MSE = c(your_mse, naive_mean_mse, random_walk_mse, ar1_mse)
)

print(mse_results) # Something is probably wrong

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
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
  
  # Generate mock factors and loadings
  factors_list <- lapply(1:T, function(t) rnorm(m))
  factors_list <- matrix(unlist(factors_list), ncol=m, byrow=T)
  loadings_list <- lapply(1:T, function(t) matrix(rnorm(p * m), ncol = m))
  residuals <- matrix(rnorm(T * p), ncol = p)
  global_factors <- matrix(rnorm(T * m), ncol = m)
  global_loadings <- matrix(rnorm(p * m), ncol = m)
  factor_covariance <- diag(m)
  
  # Compute components
  M_hat <- compute_M_hat(factors_list, global_factors, loadings_list, global_loadings, T, p)
  B_pT <- compute_B_pT(factors_list, global_factors, residuals, bandwidth, T, p, epanechnikov_kernel)
  V_pT <- compute_V_pT(factors_list, residuals, bandwidth, T, p, factor_covariance, epanechnikov_kernel)
  J_pT <- compute_J_pT(B_pT, V_pT, M_hat, T, p, bandwidth)
  
  # Print results
  print(paste("M_hat:", M_hat))
  print(paste("B_pT:", B_pT))
  print(paste("V_pT:", V_pT))
  print(paste("J_pT:", J_pT))
}


# Test forecasting implementation
test_forecasting <- function() {
  set.seed(123)
  T <- 200
  p <- 50
  W <- 100
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
  forecasts <- matrix(NA, nrow = T, ncol = p)
  
  for (tau in (W + 1):T) {
    weight_vec <- sapply(1:(tau - 1), function(r) boundary_kernel(tau - 1, r, tau - 1, bandwidth, epanechnikov_kernel))
    weight_vec <- weight_vec / sum(weight_vec)
    sqrt_w <- sqrt(weight_vec)
    weighted_subreturns <- sweep(returns[1:(tau - 1), , drop = FALSE], 1, sqrt_w, `*`)
    pca_result <- prcomp(weighted_subreturns, center = FALSE, scale = FALSE)
    local_factors_mat <- pca_result$x[, 1:m]
    loadings_hat <- t(local_factors_mat) %*% weighted_subreturns / (tau - 1)
    loadings_hat <- t(loadings_hat)
    F_tau_hat <- local_factors_mat[nrow(local_factors_mat), ]
    Rhat_tau <- loadings_hat %*% F_tau_hat
    forecasts[tau, ] <- Rhat_tau
  }
  
  # Check forecast dimensions and print a sample forecast
  print(dim(forecasts))
  print(forecasts[(W + 1):(W + 5), 1:5])  # First few forecasts for first 5 assets
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
