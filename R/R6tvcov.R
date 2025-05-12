# time_varying_cov

#time_varying_cov <- function(returns,
TVMVP$set("public", "time_varying_cov", function(
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    M0 = 10,
    rho_grid = seq(0.005, 2, length.out = 30),
    floor_value = 1e-12,
    epsilon2 = 1e-6,
    full_output = FALSE
    ) {
  flag = TRUE
  if(is.null(private$data)) {
    cli::cli_alert_warning("data is empty")
    flag = FALSE
  }
  if(!flag) return(NULL) # return

  # by default use the optimal m stored in the object, otherwise
  # determine optimal number of factors using Silverman’s bandwidth
  if(is.null(private$optimal_m)){
    m <- self$determine_factors(max_m=max_factors)
  } else {
    m <- private$optimal_m
  }
  # there should be bandwidth
  bandwidth <- private$bandwidth

  # Step 1: Local PCA
  local_res <- localPCA(
    returns     = private$data,
    bandwidth   = bandwidth,
    m           = m,
    kernel_func = kernel_func
  )

  # Step 2: Residual covariance estimation with POET
  res <- estimate_residual_cov_poet_local(
    localPCA_results = local_res,
    returns           = private$data,
    M0                = M0,
    rho_grid          = rho_grid,
    floor_value       = floor_value,
    epsilon2          = epsilon2
  )

  if (full_output) {
    return(res)
  } else {
    return(res$total_cov)
  }
})
