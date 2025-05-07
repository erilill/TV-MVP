# hypothesis testing

TVMVP$set("public", "hyptest", function(iB = 200, kernel_func = epanechnikov_kernel) {
#hyptest1 <- function(returns, m, B = 200, kernel_func = epanechnikov_kernel) {
  flag = TRUE
  if(is.null(private$data)) {
    cli::cli_alert_warning("data is empty")
    flag = FALSE
  }
  if(is.null(private$optimal_m)) {
    cli::cli_alert_warning("run {.code determine_factors()} first")
    flag = FALSE
  }
  if(!flag) return(invisible(self)) # return

  # Standardize returns
  returns <- scale(private$data)
  iT = private$iT; ip = private$ip
  if(is.null(private$bandwidth)) {
    # ??? shall we always use Silverman here?
    self$silverman()
    cli::cli_alert_warning("use Silverman bandwidth")
  }
  h = private$bandwidth
  m = private$optimal_m
  # if so then
  # h <- silverman(returns)

  # Local PCA
  localPCA_results <- localPCA(returns, bandwidth = h, m = m)
  local_factors <- localPCA_results$f_hat
  local_loadings <- localPCA_results$loadings

  # Global factor analysis
  my_svd_global <- svd(returns, nu = m, nv = m)  # Only compute top-m components
  U_m <- my_svd_global$u   # T x m
  D   <- my_svd_global$d   # length(min(T,p))
  V_m <- my_svd_global$v   # p x m

  # Match local style:
  F_global <- sqrt(iT) * U_m              # T x m
  B_global <- t((1/iT) * t(F_global) %*% returns)  # (p x m)


  # Residuals
  res <- residuals(local_factors, local_loadings, returns)
  sigma_0 <- compute_sigma_0(res, iT, ip)

  # Compute test statistic
  M_hat <- compute_M_hat(local_factors, F_global, local_loadings, B_global, iT, ip, m)
  B_pT <- compute_B_pT(local_factors, F_global, res, h, iT, ip, kernel_func)
  V_pT <- compute_V_pT(local_factors, res, h, iT, ip, kernel_func)
  J_pT <- (iT * sqrt(ip) * sqrt(h) * M_hat - B_pT) / sqrt(V_pT)

  # Step 2-4: Bootstrap procedure using sapply()
  J_pT_bootstrap <- sapply(1:iB, function(b) {
    # Step 2: Generate bootstrap error e*_it
    zeta_star <- matrix(rnorm(iT * ip, mean = 0, sd = 1), nrow = iT, ncol = ip)  # IID N(0,1)
    e_star <- t(sqrt_matrix(sigma_0) %*% t(zeta_star))  # Ensure T × p

    # Step 3: Generate new sample X*_it
    X_star <- F_global %*% t(B_global) + e_star

    # Re-run PCA for bootstrapped data
    svd_star <- svd(X_star, nu = m, nv = m)
    F_global_star <- sqrt(iT) * svd_star$u
    B_global_star <- t((1/iT) * t(F_global_star) %*% X_star)  # p × m

    # Re-run local PCA for bootstrapped data
    star_local_PCA <- localPCA(X_star, h, m)
    local_factors_star <- star_local_PCA$f_hat
    local_loadings_star <- star_local_PCA$loadings

    # Compute new residuals
    res_star <- residuals(local_factors_star, local_loadings_star, X_star)

    # Compute bootstrap test statistic J_pT*
    M_hat_star <- compute_M_hat(local_factors_star, F_global_star, local_loadings_star, B_global_star, iT, ip, m)
    B_pT_star <- compute_B_pT(local_factors_star, F_global_star, res_star, h, iT, ip, kernel_func)
    V_pT_star <- compute_V_pT(local_factors_star, res_star, h, iT, ip, kernel_func)
    return((iT * sqrt(ip) * sqrt(h) * M_hat_star - B_pT_star) / sqrt(V_pT_star))
  })
  J_pT_bootstrap <- as.numeric(unlist(J_pT_bootstrap))
  J_pT <- as.numeric(J_pT)

  # Step 4: Compute bootstrap p-value
  p_value <- mean(J_pT_bootstrap >= J_pT)

  if (p_value < 0.05) {
    message(sprintf("J_pT = %.4f, p-value = %.4f: Strong evidence that the covariance is time-varying.", J_pT, p_value))
  } else if (p_value < 0.10) {
    message(sprintf("J_pT = %.4f, p-value = %.4f: Some evidence of time-variation, but not strong.", J_pT, p_value))
  } else {
    message(sprintf("J_pT = %.4f, p-value = %.4f: No significant evidence of time-varying covariance.", J_pT, p_value))
  }

  # return(list(J_NT = J_pT, p_value = p_value, J_pT_bootstrap = J_pT_bootstrap))

  # results and return
  private$J_test = list(J_NT = J_pT, p_value = p_value, J_pT_bootstrap = J_pT_bootstrap)
  invisible(self)
})
