library(MASS)       # For matrix operations
library(matrixcalc) # For matrix calculations
library(glmnet)     # For lasso penalization
library(FactoMineR) # For PCA
library(quadprog)   # For quadratic programming
library(FE)         # Testing with real data
library(forecast)   # For ARIMA comparison

################################################################################

# Example with simulated data
set.seed(123)
T <- 200  # Number of time periods
p <- 50   # Number of assets
returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

# Example with real data (Warning: Slow)
#returns <- portfolio_m[,5:124]
#T <- nrow(returns)
#p <- ncol(returns)

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

# Function to estimate IC(m) and find the optimal number of factors
compute_V_m <- function(returns, m, kernel_func, bandwidth) {
  p <- ncol(returns)
  T <- nrow(returns)
  total_ssr <- 0  # will accumulate sum of squared residuals

  for (x in 1:T) {
    # 1) Build local weights w_x(t) around target time x
    w_x <- sapply(1:T, function(t) kernel_func((x - t)/(T * bandwidth))) 
    # or boundary_kernel(...), possibly normalized so sum(w_x)=1
    w_x <- w_x / sum(w_x)

    # 2) Form the weighted returns for local PCA
    #    WeightedRet is T x p, weighting each row t by sqrt(w_x[t])
    sqrt_w_x <- sqrt(w_x)
    WeightedRet <- sweep(returns, 1, sqrt_w_x, `*`)

    # 3) Local PCA -> factor matrix (size T x m)
    pca_result <- prcomp(WeightedRet, center = FALSE, scale = FALSE)
    num_pcs <- ncol(pca_result$x)
    m_actual <- min(m, num_pcs)
    Fhat <- pca_result$x[, 1:m_actual, drop=FALSE] # factor scores for ALL times t=1..T
    # (Optionally scale factors by 1/sqrt(T), depending on your convention)

    # 4) Weighted LS regression to get loadings B_x
    #    We treat each column of WeightedRet as an outcome, and WeightedF is design
    WeightedF  <- sweep(Fhat, 1, sqrt_w_x, `*`)        # T x m
    WeightedY  <- WeightedRet                           # T x p
    # Solve (WeightedF^T WeightedF + ridge I) B_x = WeightedF^T WeightedY
    matA <- crossprod(WeightedF)                        # (m x m)
    matB <- crossprod(WeightedF, WeightedY)             # (m x p)
    B_x  <- solve(matA, matB)                           # (m x p)

    # 5) Compute unweighted residuals and their sum of squares
    #    Resid_x(t) = returns[t, ] - (Fhat[t, ] %*% B_x)
    #    This is T x p
    fitted    <- Fhat %*% B_x                           # (T x p)
    Resid_x   <- returns - fitted
    total_ssr <- total_ssr + sum(Resid_x^2)
  }

  # 6) Finally V_m = (1 / (p*T)) * sum of squared residuals
  V_m <- total_ssr / (p * T)
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
  # Weights
  w_r <- sapply(1:T, function(t) boundary_kernel(r, t, T, bandwidth, kernel_func))
  w_r <- w_r / sum(w_r)
  
  # Weighted returns
  sqrt_w_r <- sqrt(w_r)
  weighted_returns <- sweep(returns, 1, sqrt_w_r, `*`)
  
  # Local PCA
  pca_result <- prcomp(weighted_returns, center = FALSE, scale. = FALSE)
  # Possibly check rank
  m_actual <- min(m, ncol(pca_result$x))  
  Fhat_all <- pca_result$x[, 1:m_actual, drop=FALSE]/sqrt(T)  # T x m_actual
  
  # loadings is (p x m_actual)
  loadings_all <- pca_result$rotation[, 1:m_actual, drop=FALSE]
  
  # Return the entire factor matrix
  list(
    factors_full  = Fhat_all,
    loadings_full = loadings_all
  )
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

compute_M_hat <- function(local_factors, global_factors, local_loadings, global_loadings, T, N) {
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
  B_pT <- 0
  if (m ==1){
    global_factors <- matrix(global_factors)
  }
  for (i in 1:p) {
    for (t in 1:T) {
      for (s in 1:T) {
        k_h_st <- boundary_kernel(s, t, T, h, kernel_func)
        diff <- (k_h_st * (t(local_factors[[s]][s,]) %*% local_factors[[t]][t,]) - # Outer poroduct or dot product?
                   (t(global_factors[s, ]) %*% global_factors[t, ]))^2
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
M_hat <- compute_M_hat(factors_list, global_factors, loadings_list, global_loadings, T, p)
print(paste("M_hat:", M_hat))

# Compute B_pT
B_pT <- compute_B_pT(factors_list, global_factors, residuals, bandwidth, T, p, epanechnikov_kernel)
print(paste("B_pT:", B_pT))

# Compute V_pT
V_pT <- compute_V_pT(factors_list, residuals, bandwidth, T, p, factor_covariance, epanechnikov_kernel)
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
    boundary_kernel(tau - 1, r, (tau - 1), bandwidth, epanechnikov_kernel)
  })
  weight_vec <- weight_vec / sum(weight_vec)  # Normalize weights
  
  # Weighted returns: up to time tau-1
  sqrt_w <- sqrt(weight_vec)
  weighted_subreturns <- sweep(returns[1:(tau-1), , drop = FALSE], 1, sqrt_w, `*`)
  
  # Perform PCA
  pca_result <- prcomp(weighted_subreturns, center = FALSE, scale = FALSE)
  
  # Check the number of components returned by PCA
  num_factors <- min(m, ncol(pca_result$x))  # Ensure we don't exceed available factors
  if (num_factors == 0) {
    warning(paste("Skipping tau =", tau, "due to insufficient factors."))
    next
  }
  
  local_factors_mat <- pca_result$x[, 1:num_factors, drop = FALSE]  # (tau-1) x num_factors
  loadings_hat <- t(local_factors_mat) %*% weighted_subreturns / (tau - 1)
  loadings_hat <- t(loadings_hat)  # p x num_factors
  
  # Handle the case where num_factors = 1 (ensure factors are treated as a matrix)
  if (num_factors == 1) {
    F_tau_hat <- local_factors_mat[nrow(local_factors_mat), drop = FALSE]  # Ensure it's a matrix
  } else {
    F_tau_hat <- local_factors_mat[nrow(local_factors_mat), ]  # Last row
  }
  
  # Forecast
  Rhat_tau <- loadings_hat %*% F_tau_hat
  
  # Store forecast for time tau
  forecasts[tau, ] <- Rhat_tau
}


# Initialize matrices for alternative forecasts
naive_mean_forecasts <- matrix(NA, nrow = T, ncol = p)
random_walk_forecasts <- matrix(NA, nrow = T, ncol = p)
ar1_forecasts <- matrix(NA, nrow = T, ncol = p)

# Loop over assets
fit_ar1 <- ar(returns[1:W, p], order.max = 1, method = "ols")

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
      ar1_forecasts[tau, j] <- predict(fit_ar1, newdata = returns[1:tau, j], n.ahead = 1)$pred
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

