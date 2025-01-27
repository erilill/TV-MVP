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
# - local_pca, determine_factors, localPCA: Follows Su & Wang now but needs to be optimized
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
m <- determine_factors(returns, ncol(returns)) # Needs optimization
local_pca_res <- localPCA(returns, bandwidth, m$optimal_R) # Needs optimization
summary(local_pca_res$factors)
t(local_pca_res$factors)%*%local_pca_res$factors*(1/nrow(returns)) #covariance

# Global PCA
global_pca <- prcomp(returns, scale. = FALSE, center = TRUE)
global_factors <- global_pca$x[, 1:m$optimal_R]
global_loadings <- global_pca$rotation[, 1:m$optimal_R]

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
