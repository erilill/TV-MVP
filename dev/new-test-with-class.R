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

tmp$set(data = returns)
tmp$set_data(returns)
tmp
# show head like tibble does

#tmp$set_data(tibble::as_tibble(returns))
#tmp$get_data()

# this does not work because x y a b are not defined in the class
#tmp$set(x=1, y=2)$set(a=3)$set(b=4)
# but this will work
#tmp$set(iT=10)$set(ip=2)

# warning when first run
tmp$determine_factors()
# works, use the default Silverman
tmp$determine_factors(10)
tmp
tmp$get_optimal_m()
tmp$get_IC_values()
# I think the output should be simply optimal_m, then perhaps the user could use:
# tmp$IC in order to get the IC-values in case they want to see them or plot them.
## done! the function returns optimal_m,
## user can get the optimal_m by using get_optimal_m
## we don't do active binding...

tmp$hyptest(iB = 10, kernel_func = epanechnikov_kernel)
tmp
tmp$get_bootstrap()
# Here, the output should be the test statistic (first value in the list),
# and the bootstrap p-value (second value in the list).
# The vector contains all of the bootstrap test statistics which might be of
# interest to the user if they want to plot it, but it is likely not necessary
# to show at first glance. Similarly, we could include something like:
# tmp$J_star which prints these.

# some bug here... does not work
prediction <- predict_portfolio(returns, horizon = 21, epanechnikov_kernel, m = 10, min_return=0.5)
# it works here, weird...
prediction <- tmp$predict_portfolio(horizon = 1,
                                    kernel_func = epanechnikov_kernel,
                                    min_return = 0.01,
                                    max_SR = TRUE,
                                    rf = 0.001)
prediction
prediction$getWeights(method = "MVP")
prediction$getWeights(method = "max_SR")
prediction$getWeights(method = "MVPConstrained")
# This function could probably quite easily be included in the class if we remove
# the computation of max_m within the function or set a condition so that it does
# not run determine_factors() if the class already has max_m computed.
# I think this would be nice to include in the class, perhaps with default settings:
# bandwidth = silverman, min_return = NULL, max_SR = TRUE. The user then runs
# tmp$predict(horizon = 21) or tmp$predict(horizon = 21, min_return = 0.5).

# Evaluate historical performance of model:
mvp_result <- rolling_time_varying_mvp(
  returns        = returns,
  initial_window = 60,
  rebal_period   = 5,
  max_factors    = 10,
  return_type    = "daily",
  kernel_func    = epanechnikov_kernel
)
mvp_result

mvp_result <- tmp$rolling_time_varying_mvp(
  initial_window = 60,
  rebal_period   = 5,
  max_factors    = 10,
  return_type    = "daily",
  kernel_func    = epanechnikov_kernel
)
mvp_result
# I like this function and would like to include it in the package, however I
# am not certain of how it fits in the class. There are quite a lot of parameters
# to be set which cant really have a default: initial_window depends on size of
# dataset, rebal_period depends on the aggregation level of the data and intended
# investment horizon, and the return_type also depends on the aggregation level.

# Further, this function cannot use the already computed optimal_m as it computes
# optimal_m based on the initial window.

# I think rolling_time_varying_mvp is the most interesting to plot, as this will
# show the portfolio performance over time.


# I also think time_varying_cov() should be included in the class:
cov_mat <- time_varying_cov(returns,
                                m=1,
                                bandwidth = silverman(returns),
                                kernel_func = epanechnikov_kernel,
                                M0 = 10,
                                rho_grid = seq (0.005 , 2,
                                                  length.out = 30),
                                floor_value = 1e-12,
                                epsilon2 = 1e-6,
                                full_output = FALSE)
cov_mat

cov_mat <- tmp$time_varying_cov(kernel_func = epanechnikov_kernel,
                                M0 = 10,
                                rho_grid = seq (0.005 , 2,
                                                length.out = 30),
                                floor_value = 1e-12,
                                epsilon2 = 1e-6,
                                full_output = FALSE)
cov_mat
# full_output = TRUE prints some additional information, this could be of interest
# to some users who wants more information. However, the most important output is
# the covariance matrix, "total_cov". The user could then construct their portfolio
# using this.
