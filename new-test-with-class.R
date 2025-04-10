# Load necessary libraries
library(TVMVP)

# Things that I am unsure of/needs work:
# - determine_factors, localPCA: slow
# - cv_bandwidth: is this necessary? Might need to look at other articles.
# - rolling_time_varying_mvp, predict_portfolio: I think it works
# - compute_residual_covariance: Works but problems with singularity, probably due to some problem with local_pca
# - local_pca: Worked a lot in order to get this to work properly, still not sure.
# - Licence of quadprog and spcov: GPL 2/>=2. If I want to use them, do I also need to use GNU?
# - Realised that it is not really necessary for the end user to to do the localPCA themselves. Made it more userfriendly.

set.seed(123)
T <- 100  # Number of time periods
p <- 20   # Number of assets
returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

tmp <- TVMVP$new()
tmp$set_data()
tmp$set_data(tibble::as_tibble(returns))
tmp$get_data()

# Number of factors
m <- determine_factors(returns, 10, silverman(returns))$optimal_R # Needs optimization

# Test if covariance is time invariant
hyptest1(returns = returns,
         m,
         B = 199,
         kernel_func = epanechnikov_kernel
)
# E: Slow, but I think it works. The test statistics are sometimes much larger than expected. But p-vals seems correct.

# Evaluate historical performance of model:
mvp_result <- rolling_time_varying_mvp(
  returns        = returns,
  initial_window = 60,
  rebal_period   = 5,
  max_factors    = 10,
  return_type    = "daily",
  kernel_func    = epanechnikov_kernel,
  bandwidth_func = silverman
)
## K: some problem with non-conformable arguments
## K: Error in factors_t %*% t(loadings_t) : non-conformable arguments
## E: When initial_window is small, and m is large, the effective rank in the cv_bandwidth becomes < m which causes problems
## E: Tried to fix it. I am not sure if this is a valid solution to the problem.

prediction <- predict_portfolio(returns, horizon = 21, silverman, max_factors = 10, min_return=0.5)
