# Load necessary libraries
library(TVMVP)

# Things that I am unsure of/needs work:
# - determine_factors, localPCA: slow
# - rolling_time_varying_mvp, predict_portfolio: I think it works
# - local_pca: Worked a lot in order to get this to work properly, still not sure.

set.seed(123)
uT <- 100  # Number of time periods
up <- 20   # Number of assets
returns <- matrix(rnorm(uT * up, mean = 0.001, sd = 0.02), ncol = up)

# Number of factors
m <- determine_factors(returns, 10, silverman(returns))$optimal_m # Needs optimization

# Test if covariance is time invariant
test <- hyptest1(returns = returns,
         m,
         B = 199,
         kernel_func = epanechnikov_kernel)
# E: Slow, but I think it works. The test statistics are sometimes much larger than expected. But p-vals seems correct.

# Evaluate historical performance of model:
mvp_result <- rolling_time_varying_mvp(
  returns        = returns,
  initial_window = 60,
  rebal_period   = 5,
  max_factors    = 10,
  return_type    = "daily",
  kernel_func    = epanechnikov_kernel,
  rf             = 1e-04
)
prediction <- predict_portfolio(returns, 
                                horizon = 21, 
                                max_factors = 10,
                                kernel_func = epanechnikov_kernel,
                                min_return=0.5,
                                max_SR = TRUE,
                                rf = 1e-04)
