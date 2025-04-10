## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(TVMVP)

## -----------------------------------------------------------------------------
set.seed(123)
uT <- 100  # Number of time periods
up <- 20   # Number of assets
returns <- matrix(rnorm(uT * up, mean = 0.001, sd = 0.02), ncol = up)

## -----------------------------------------------------------------------------
m <- determine_factors(returns = returns, max_m = 10, bandwidth = silverman(returns))$optimal_m
m
hypothesis_test <- hyptest1(returns = returns,
                            m = m,
                            B = 10, # Use larger B in practice
                            kernel_func = epanechnikov_kernel)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
prediction <- predict_portfolio(returns = returns, 
                                horizon = 21, 
                                max_factors = 10,
                                kernel_func = epanechnikov_kernel,
                                min_return=0.5,
                                max_SR = TRUE,
                                rf = 1e-04)
prediction

## -----------------------------------------------------------------------------
cov_mat <- time_varying_cov(returns,
                            m,
                            bandwidth = silverman(returns),
                            kernel_func = epanechnikov_kernel,
                            M0 = 10,
                            rho_grid = seq(0.005, 2, length.out = 30),
                            floor_value = 1e-12,
                            epsilon2 = 1e-6,
                            full_output = FALSE)
  

