# Load necessary libraries
library(TVMVP)

# Things that I am unsure of/needs work:
# - determine_factors, localPCA: slow
# - rolling_time_varying_mvp, predict_portfolio: I think it works
# - local_pca: Worked a lot in order to get this to work properly, still not sure.

set.seed(123)
simulate_time_varying_data <- function(T, p, phi, sigma_f){
  F       <- numeric(T)
  F[1] <- rnorm(1, mean = 0, sd = sigma_f / sqrt(1 - phi^2))
  for (t in 2:T) {
    sigma_ft <- sigma_f
    F[t] <- phi * F[t - 1] + rnorm(1, 0, sigma_ft)
    }
  alpha    <- runif(p, min = 0.0002, max = 0.001)
  loadings <- runif(p, min = 0.5,    max = 1.5)
  returns <- matrix(NA, nrow = T, ncol = p)
  for (i in 1:p) {
    for (t in 1:T) {
      sigma_e_t <- 0.01 + 0.005 * sin(2 * pi * t / 50 + i * 0.1)
      e_it <- rnorm(1, mean = 0, sd = sigma_e_t)
      returns[t, i] <- alpha[i] + loadings[i] * F[t] + e_it
    }
  }
  return(returns)
}
returns <- simulate_time_varying_data(100, 30, 0.7, 0.01)

# Number of factors
m <- determine_factors(returns, 10, silverman(returns))$optimal_m # Needs optimization
m

# Test if covariance is time invariant
test <- hyptest1(returns = returns,
         m,
         B = 199,
         kernel_func = epanechnikov_kernel)
test
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
mvp_result
plot(mvp_result)
prediction <- predict_portfolio(returns, 
                                horizon = 21, 
                                max_factors = 10,
                                kernel_func = epanechnikov_kernel,
                                min_return=0.5,
                                max_SR = TRUE,
                                rf = 1e-04)
prediction