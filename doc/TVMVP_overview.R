## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----example------------------------------------------------------------------
library(TVMVP)

## ----data---------------------------------------------------------------------
set.seed(123)
uT <- 100  # Number of time periods
up <- 20   # Number of assets
returns <- matrix(rnorm(uT * up, mean = 0.001, sd = 0.02), ncol = up)

## ----initialize---------------------------------------------------------------
tvmvp_obj <- TVMVP$new()
tvmvp_obj$set_data(returns)

## ----hyptest------------------------------------------------------------------
tvmvp_obj$determine_factors(max_m=5)
tvmvp_obj$get_optimal_m()

tvmvp_obj$hyptest(iB=10) # Use larger iB in practice
tvmvp_obj

## ----rolpred------------------------------------------------------------------
mvp_result <- tvmvp_obj$rolling_time_varying_mvp(
  initial_window = 60,
  rebal_period   = 5,
  max_factors    = 10,
  return_type    = "daily",
  rf             = NULL
)

mvp_result

## ----pred---------------------------------------------------------------------
prediction <- tvmvp_obj$predict_portfolio(horizon = 21, min_return = 0.5, 
                                   max_SR = TRUE)
prediction
weights <- prediction$getWeights("MVP")

## ----cov----------------------------------------------------------------------
cov_mat <- tvmvp_obj$time_varying_cov()

## ----functionex, eval=FALSE---------------------------------------------------
# # Determine number of factors
# m <- determine_factors(returns = returns, max_m = 10, bandwidth = silverman(returns))$optimal_m
# m
# 
# # Run test of constant loadings
# hypothesis_test <- hyptest1(returns = returns,
#                             m = m,
#                             B = 10, # Use larger B in practice
#                             )
# 
# # Rolling window evaluation
# mvp_result <- rolling_time_varying_mvp(
#   returns        = returns,
#   initial_window = 60,
#   rebal_period   = 5,
#   max_factors    = 10,
#   return_type    = "daily",
#   kernel_func    = epanechnikov_kernel,
#   rf             = 1e-04
# )
# mvp_result
# 
# # Optimize weights and predict performance out-of-sample
# prediction <- predict_portfolio(returns = returns,
#                                 horizon = 21,
#                                 m = 10,
#                                 kernel_func = epanechnikov_kernel,
#                                 min_return=0.5,
#                                 max_SR = TRUE,
#                                 rf = 1e-04)
# prediction
# weights <- prediction$getWeights("MVP")
# 
# # For custom portfolio optimization, compute the time dependent covariance:
# cov_mat <- time_varying_cov(returns,
#                             m,
#                             bandwidth = silverman(returns),
#                             kernel_func = epanechnikov_kernel,
#                             M0 = 10,
#                             rho_grid = seq(0.005, 2, length.out = 30),
#                             floor_value = 1e-12,
#                             epsilon2 = 1e-6,
#                             full_output = FALSE)
# 

