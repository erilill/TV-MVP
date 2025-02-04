library(devtools)
library(roxygen2)
build()
document()
install()

# Load necessary libraries
library(TVMVP)

# Things that I am unsure of/needs work:
# - determine_factors, localPCA: slow
# - cv_bandwidth: is this necessary? Might need to look at other articles.
# - rolling_time_varying_mvp, predict_portfolio: I think it works
# - compute_residual_covariance: Works but problems with singularity, probably due to some problem with local_pca
# - local_pca: Worked a lot in order to get this to work properly, still not sure.
# - Fan et al. uses CV to compute lambda for the residual covariance, will look into this.

set.seed(123)
T <- 100  # Number of time periods
p <- 20   # Number of assets
returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

# Number of factors
m <- determine_factors(returns, 10, silverman(returns))$optimal_R # Needs optimization

# Select bandwidth using Silverman's rule of thumb or CV
bandwidth <- silverman(returns)

bandwidth <- cv_bandwidth(returns, m = m, candidate_h = seq(0.05, 0.95, 0.05),
                          kernel_func = epanechnikov_kernel)$optimal_h # Needs optimization
## K: Warnings about many singularities

# Perform Local PCA for all time points
local_pca_res <- localPCA(returns, bandwidth, m) # Needs optimization
summary(local_pca_res$factors)
t(local_pca_res$factors)%*%local_pca_res$factors*(1/nrow(returns)) #covariance

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
  rebal_period   = 5,
  max_factors    = 10,
  kernel_func    = epanechnikov_kernel,
  bandwidth_func = cv_bandwidth
)
## K: some problem with non-conformable arguments
## K: Error in factors_t %*% t(loadings_t) : non-conformable arguments
## E: When initial_window is small, and m is large, the effective rank in the cv_bandwidth becomes < m which causes problems
## E: Tried to fix it. I am not sure if this is a valid solution to the problem.

prediction <- predict_portfolio(returns, horizon = 21, silverman, max_factors = 10, min_return=0.5)
