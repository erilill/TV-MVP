#' @export
estimate_residual_cov_poet_local <- function(localPCA_results, 
                                             returns,
                                             M0 = 10, 
                                             rho_grid = seq(0.005, 2, length.out = 30),
                                             floor_value = 1e-12,
                                             epsilon2 = 1e-6) {

  # This function:
  #   1. Form local residuals u_t = R_local - F(z_t) * B(z_t)ᵀ
  #   2. Call adaptive_poet_rho() on those residuals to pick the single best rho_t
  #   3. Compute raw residual covariance, shrink once using rho_t
  #   4. Combine with factor part to get Sigma_X(z_t)
  #
  # The function returns a list with the final local residual and total covariances.
  

  
    # 1. Extract local loadings, factors, and row indices
    Lambda_t <- localPCA_results$loadings[[nrow(returns)]]  # p x K
    F_t <- localPCA_results$f_hat   # T_t x K
    w_t <- localPCA_results$weights[[nrow(returns)]]
    idx <- t  # subset of rows in 'returns' for this window
    
    # 2. residuals
    U_local <- residuals(F_t, localPCA_results$loadings, returns)  # T_t x p
    
    # 3. Pick best rho for these local residuals using Chen–Leng–style grouping
    #    (See the 'adaptive_poet_rho()' function you already have.)
    #    This returns (best_rho = ..., min_Fnorm = ...)
    rho_result <- adaptive_poet_rho(U_local,
                                           M0 = M0,
                                           rho_grid = rho_grid,
                                           epsilon2 = epsilon2)
    best_rho_t <- rho_result$best_rho
    
    # 4. Compute the naive residual covariance, then shrink once
    S_u_raw <- (1 / nrow(U_local)) * crossprod(U_local)  # p x p
    # threshold based on best_rho_t
    threshold <- best_rho_t * mean(abs(S_u_raw[upper.tri(S_u_raw)]))
    
    # soft-threshold function
    soft_threshold <- function(x, thr) {
      sign(x) * pmax(abs(x) - thr, 0)
    }
    # apply to off-diagonal entries
    S_u_shrunk <- apply(S_u_raw, c(1, 2), soft_threshold, thr = threshold)
    # keep diagonal as-is
    diag(S_u_shrunk) <- diag(S_u_raw)
    
    # 5. Final local covariance = factor part + shrunk residual
    Sigma_R_t <- Lambda_t%*%(crossprod(F_t)/nrow(returns))%*%t(Lambda_t) + S_u_shrunk  # p x p
    
    #Please check /Erik
    ############################################################################
    # Final PSD repair
    e_decomp <- eigen(Sigma_R_t, symmetric = TRUE)
    eigvals  <- e_decomp$values
    eigvecs  <- e_decomp$vectors
    eigvals_floored <- ifelse(eigvals < floor_value, floor_value, eigvals)
    Sigma_R_t <- eigvecs %*% diag(eigvals_floored) %*% t(eigvecs) # reconstruct sigma
    Sigma_R_t <- 0.5 * (Sigma_R_t + t(Sigma_R_t)) #Symmetrize
    ############################################################################
    
    
    # store results
    output_list <- list(
      best_rho        = best_rho_t,
      residual_cov    = S_u_shrunk,  # Σ̂e(T)
      total_cov       = Sigma_R_t,   # Σ̂_R(T)
      loadings        = Lambda_t,
      naive_resid_cov = S_u_raw 
    )
  
  
  return(output_list)
}



#' @export
adaptive_poet_rho <- function(R, M0 = 10,
                                     rho_grid = seq(0.001, 2, length.out = 20),
                                     epsilon2 = 1e-6) {
  # R: data matrix, dimension T x p
  # M0: number of observations to leave out between the two sub-samples
  # rho_grid: grid of possible rho values
  
  iT <- nrow(R)
  ip <- ncol(R)
  
  # Half of the sample size (floored)
  halfT <- floor(iT / 2)
  
  # Define T1 and T2 for the sub-samples
  T1 <- floor(halfT * (1 - 1 / log(iT)))
  T2 <- halfT - T1  # ensures T1 + T2 = floor(T/2)
  
  # Number of groups as per Chen et al. (2019), p. 61:
  # "we divide the full sample into floor(T/(2*M0)) groups"
  num_groups <- floor(iT / (2 * M0))
  
  if (num_groups < 1) {
    stop("Not enough data for adaptive rho selection. Increase T or reduce M0.")
  }
  
  # A small helper for soft-thresholding
  soft_threshold <- function(x, thr) {
    sign(x) * pmax(abs(x) - thr, 0)
  }
  
  # Function to compute the SUM of Frobenius-norm differences
  # across all groups for a given rho
  frob_sum_for_rho <- function(rho) {
    total_error <- 0
    lambda_min_vals <- numeric(num_groups)
    
    for (m in seq_len(num_groups)) {
      # The m-th group includes observations from:
      #    start_idx = (m - 1)*M0 + 1
      #    end_idx   = (m - 1)*M0 + (halfT + M0)
      #
      # each group is of length halfT + M0, with M0 overlap/step.
      
      start_idx <- (m - 1) * M0 + 1
      end_idx   <- (m - 1) * M0 + (halfT + M0)
      if (end_idx > iT) break  # guard for edge cases
      
      # Sub-sample 1: first T1 observations of the group
      sub1_start <- start_idx
      sub1_end   <- start_idx + T1 - 1
      
      # Skip M0 observations after sub1
      # sub-sample 2: last T2 observations
      sub2_start <- start_idx + T1 + M0
      sub2_end   <- sub2_start + T2 - 1
      
      if (sub2_end > end_idx) break  # guard for edge cases
      
      # Extract the actual data for sub-samples
      data_sub1 <- R[sub1_start:sub1_end, , drop = FALSE]
      data_sub2 <- R[sub2_start:sub2_end, , drop = FALSE]
      
      # Covariance from the first sub-sample: will be shrunk
      S1 <- cov(data_sub1)
      # Covariance from the second sub-sample: "naive" benchmark
      S2 <- cov(data_sub2)
      
      # Apply soft-thresholding to S1 based on rho
      threshold <- rho * mean(abs(S1[upper.tri(S1)]))
      S1_shrunk <- apply(S1, c(1, 2), soft_threshold, thr = threshold)
      
      eigvals <- eigen(S1_shrunk, symmetric = TRUE, only.values = TRUE)$values
      lambda_min_vals[m] <- min(eigvals)
      
      # Accumulate Frobenius norm difference
      total_error <- total_error + sum((S1_shrunk - S2)^2)
    }
    
    return(list(total_error = total_error, lambda_min_vals = lambda_min_vals))
  }
  
  # We now scan across rho_grid, compute the total Frobenius difference,
  # and pick the single rho that MINIMIZES the sum across all groups.
  
  best_rho <- NA
  min_val  <- Inf
  lambda_min_all <- numeric(length(rho_grid))
  for (i in seq_along(rho_grid)){
    rho_result <- frob_sum_for_rho(rho_grid[i])
    lambda_min_all[i] <- min(rho_result$lambda_min_vals, na.rm=iT)
  }
  
  # Compute rho_1
  valid_rho_indices <- which(!is.na(lambda_min_all) & lambda_min_all > 0)
  rho_1 <- if (length(valid_rho_indices) > 0) {
    epsilon2 + min(rho_grid[valid_rho_indices])
  } else {
    epsilon2  # Default to epsilon2 if no valid rho is found
  }
  
  # Compute rho
  for (rho in rho_grid[rho_grid >= rho_1]) {
    rho_result <- frob_sum_for_rho(rho)
    val <- rho_result$total_error
    
    if (val < min_val) {
      min_val  <- val
      best_rho <- rho
    }
  }
  #cat(sprintf("Chosen rho = %.5f with rho_1 = %.5f, total F-norm difference = %.5f\n", best_rho, rho_1, min_val))
  
  return(list(best_rho = best_rho, rho_1 = rho_1, min_Fnorm = min_val))
}
