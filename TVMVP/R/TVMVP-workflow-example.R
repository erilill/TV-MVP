# Example

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

data("edhec")
head(edhec)
# edhec is a time series (xts) of hedge fund indexes
returns<- as.matrix(edhec)



# Select bandwidth using Silverman's rule of thumb
bandwidth <- silverman(returns)

# Perform Local PCA for all time points
local_pca_res <- localPCA(returns)

m <- 1

# Global PCA
global_pca <- prcomp(returns, scale. = FALSE, center = FALSE)
global_factors <- global_pca$x[, 1:m]/sqrt(T)
global_loadings <- global_pca$rotation[, 1:m]

# Compute residuals
res <- residuals(local_pca_res$factors_list, local_pca_res$loadings_list, returns)

# Test if covariance is time invariant
hyptest1(
  local_factors = local_pca_res$factors_list,
  global_factors = global_factors,
  local_loadings = local_pca_res$loadings_list,
  global_loadings = global_loadings,
  residuals = res,
  kernel_func = epanechnikov_kernel
)

W <- 100

# Evaluate model on data to see if it is sensible to use:
forecast <- forecast_local_factor_model(returns, W, m, bandwidth)

covs <- realized_covariances(W, returns, res)
true_weights <- optimal_weights(W, local_pca_res$factors_list,
                                local_pca_res$loadings_list,
                                covs$realized_resid_covariances,
                                nrow(returns),
                                ncol(returns))
sharpes <- SR(W, returns, local_pca_res$loadings_list,
              local_pca_res$factors_list, true_weights,
              covs$realized_resid_covariances)

ev <- evaluate_forecasts(
  returns        = returns,
  forecasts      = forecast$forecasts,
  est_covariances = forecast$est_covariances,
  residual_covariances_pred = forecast$residual_covariances,
  weights_est    = forecast$optimal_weights,
  window_eval    = (W+1):T,
  realized_covariances       = covs$realized_covariances,
  realized_resid_covariances = covs$realized_resid_covariances,
  true_weights   = true_weights,
  realized_sharpes = sharpes
)
print(ev)

# If satisfactory, use this method in order to produce the minimum variance portfolio
tvmvp <- minvar_portfolio(colMeans(returns),
                          compute_time_varying_cov(local_pca_res$loadings_list[[200]],
                                                   cov(local_pca_res$factors_list[[200]]),
                                                   covs$realized_resid_covariances[[200]]
                                                   ),
                          0.005)
plot_portfolio_weights(tvmvp$w_g, title = "Global Minimum Portfolio Weights")
plot_portfolio_weights(tvmvp$w_p, title = "Target Minimum Portfolio Weights")

