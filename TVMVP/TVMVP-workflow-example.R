library(devtools)
remove.packages("TVMVP")
devtools::clean_dll()
build()
install()

# Example
setwd("C:/Users/erikl_xzy542i/Documents/Master_local/Thesis/TV-MVP/TVMVP")
# Load necessary libraries
library(MASS)        # For matrix operations
library(matrixcalc)  # For matrix calculations
library(glmnet)      # For lasso penalization
library(FactoMineR)  # For PCA
library(quadprog)    # For quadratic programming
library(forecast)    # For ARIMA comparison
library(zoo)         # For rolling operations
library(ggplot2)     # For visualization
library(TVMVP)
library(PerformanceAnalytics)

# Things that I am unsure of/needs work:
# - localPCA: unsure if correctly computed
# - estimate_residual_covariance: Wrong, do not understand the article
# - cv_bandwidth: most likely wrong
# - mvp_result: The results are strange, could be wrong, or could be affected by other bad functions
# - predict_portfolio: Results are a bit strange, have only gotten it to work with global min var portf.



data("edhec")
head(edhec)
returns<- as.matrix(edhec)



# Select bandwidth using Silverman's rule of thumb or CV
bandwidth <- silverman(returns)

folds <- list(
  returns[1:58, ],    # Fold 1
  returns[59:118, ],   # Fold 2
  returns[119:176, ],  # Fold 3
  returns[177:236, ], # Fold 4
  returns[237:293, ]  # Fold 5
)
bandwidth <- cv_bandwidth(returns, folds, seq(0.1, 1.0, by=0.05), 10, epanechnikov_kernel)

# Perform Local PCA for all time points
local_pca_res <- localPCA(returns, bandwidth, 10)
summary(local_pca_res$factors)
cov(local_pca_res$factors)
m <- local_pca_res$m

# Global PCA
global_pca <- prcomp(returns, scale. = FALSE, center = TRUE)
global_factors <- global_pca$x[, 1:m]
global_loadings <- global_pca$rotation[, 1:m]

# Compute residuals
res <- residuals(local_pca_res$factors, local_pca_res$loadings, returns)

# Test if covariance is time invariant
hyptest1(
  local_factors = local_pca_res$factors,
  global_factors = global_factors,
  local_loadings = local_pca_res$loadings,
  global_loadings = global_loadings,
  residuals = res,
  kernel_func = epanechnikov_kernel
)

# Evaluate historical performance of model:
mvp_result <- rolling_time_varying_mvp(
  returns        = returns,
  initial_window = 60,
  rebal_period   = 20,
  max_factors    = 3,
  bandwidth      = 0.2,
  kernel_func    = epanechnikov_kernel
)

prediction <- predict_portfolio(returns, horizon = 21, bandwidth = 0.2, max_factors = 3)
























determine_factors <- function(returns, max_R, g_func = function(N, T) log(N * T) / (N * T)) {
  T <- nrow(returns)
  N <- ncol(returns)

  # Initialize storage
  V <- numeric(max_R)
  penalty <- numeric(max_R)

  # Loop over possible number of factors (R)
  for (R in 1:max_R) {
    residuals <- matrix(NA, nrow = T, ncol = N)
    for (r in 1:T){
    # Step 1: Perform PCA with R factors
    pca_result <- local_pca(returns, r = r, bandwidth = bandwidth, m = R, kernel_func = epanechnikov_kernel)
    X_r <- matrix(0, nrow = T, ncol = N)
    X_r <- sweep(returns, 1, sqrt(pca_result$w_r), `*`)
    Lambda_breve_R <- t((1/T*N)*t(X_r)%*%X_r%*%pca_result$loadings)
    F_breve_R <- solve((Lambda_breve_R)%*%t(Lambda_breve_R))%*%(Lambda_breve_R)%*%returns[r,]

    # Step 2: Compute SSR (Sum of Squared Residuals)
    residuals[r,] <- returns[r,] - t(F_breve_R) %*% (Lambda_breve_R)
    V[R] <- sum(residuals^2) / (N * T)

    penalty[R] <- R * g_func(N, T)
  }
  }
  # Step 4: Determine optimal number of factors
  IC_values <- log(V) + penalty
  optimal_R <- which.min(IC_values)

  return(list(optimal_R = optimal_R, IC_values = IC_values))
}


localPCA <- function(returns,
                     bandwidth,
                     max_factors,
                     kernel_func = epanechnikov_kernel) {
  p <- ncol(returns)
  T <- nrow(returns)

  # Example: user-supplied function to pick the number of factors
  # (replace with your actual factor selection code)
  m_selection <- select_optimal_factors(
    returns = returns,
    max_factors = max_factors,
    T_h = T * bandwidth,
    kernel_func = kernel_func,
    bandwidth = bandwidth
  )
  m <- m_selection$optimal_m

  # Initialize storage
  factors <- matrix(NA, nrow = T, ncol = m)
  loadings <- vector("list", T)
  weights_list <- vector("list", T)

  # For each time t, do local PCA
  for (t_i in 1:T) {
    local_result <- local_pca(returns, t_i, bandwidth, m, kernel_func)
      factors[t_i, ] <- local_result$factors
      loadings[[t_i]] <- local_result$loadings
      weights_list[[t_i]] <- local_result$w_r
  }

  return(list(
    factors = factors,    # T x m
    loadings = loadings,  # list of length T, each p x m
    m = m,
    weights = weights_list
  ))
}
