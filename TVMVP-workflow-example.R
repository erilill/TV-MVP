# Load necessary libraries
library(TVMVP)

# Things that I am unsure of/needs work:
# - determine_factors, localPCA: slow
# - cv_bandwidth: is this necessary? Might need to look at other articles.
# - rolling_time_varying_mvp, predict_portfolio: I think it works
# - compute_residual_covariance: Works but problems with singularity, probably due to some problem with local_pca
# - local_pca: Worked a lot in order to get this to work properly, still not sure.
# - Fan et al. uses CV to compute lambda for the residual covariance, will look into this.
# - Licence of quadprog and spcov: GPL 2/>=2. If I want to use them, do I also need to use GNU?

set.seed(123)
T <- 100  # Number of time periods
p <- 20   # Number of assets
returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

# Number of factors
m <- determine_factors(returns, 10, silverman(returns))$optimal_R # Needs optimization

# Select bandwidth using Silverman's rule of thumb or CV
bandwidth <- silverman(NULL, T, p)

bandwidth <- cv_bandwidth(returns, m = m, candidate_h = seq(0.05, 0.95, 0.05),
                          kernel_func = epanechnikov_kernel)$optimal_h # Needs optimization
## K: Warnings about many singularities
## E: Should be fixed now. I believe it was because spcov was not in the imports.
## E: Have not stress-tested it fully, but it is not stable when T is small or m is large.

# Perform Local PCA for all time points
local_pca_res <- localPCA(returns, bandwidth, m) # Needs optimization

# Test if covariance is time invariant
hyptest1(localPCA_results = local_pca_res,
         returns = returns,
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
