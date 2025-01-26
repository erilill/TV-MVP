#' Work in progress
#'
compute_optimal_weights <- function(returns_train) {
  cov_mat <- cov(returns_train)
  inv_cov <- solve(cov_mat)
  ones <- rep(1, ncol(returns_train))
  w_unnorm <- inv_cov %*% ones
  w <- as.numeric(w_unnorm / sum(w_unnorm))

  return(w)
}

realized_sharpe <- function(w, test_data) {
  port_ret <- test_data %*% w
  sr <- mean(port_ret) / sd(port_ret)
  return(sr)
}

compute_residuals <- function(local_factors, loadings_list, returns) {
  Tn <- nrow(returns)
  p <- ncol(returns)
  res <- matrix(NA, Tn, p)
  for (tt in 1:Tn) {
    modeled <- local_factors[tt, , drop=FALSE] %*% t(loadings_list[[tt]])
    res[tt, ] <- returns[tt, ] - modeled
  }
  res
}


#' @export
cv_bandwidth <- function(returns, folds, candidate_h, max_factors, kernel_func) {
  k <- length(folds)
  scores <- numeric(length(candidate_h))

  for (h_i in seq_along(candidate_h)) {
    h_val <- candidate_h[h_i]
    sr_sum <- 0
    for (j in 2:k) {
      train_data <- do.call(rbind, folds[1:(j-1)])
      test_data  <- folds[[j]]
      local_pca_train <- localPCA(train_data, bandwidth = h_val,
                                  max_factors = max_factors, kernel_func = kernel_func)
      factor_cov <- cov(local_pca_train$factors)
      res <- compute_residuals(
        local_factors = local_pca_train$factors,
        loadings_list = local_pca_train$loadings,
        returns       = train_data
      )
      residual_cov <- cov(res)
      load_sum <- matrix(0, nrow=p, ncol=m)
      valid_count <- 0
      for (tt in 1:nrow(train_data)) {
        if (!anyNA(local_pca_train$loadings[[tt]])) {
          load_sum <- load_sum + local_pca_train$loadings[[tt]]
          valid_count <- valid_count + 1
        }
      }
      avg_loadings <- load_sum / valid_count

      Sigma_hat <- avg_loadings %*% factor_cov %*% t(avg_loadings) + residual_cov
      inv_cov <- solve(Sigma_hat)
      ones <- rep(1, nrow(Sigma_hat))
      w_unnorm <- inv_cov %*% ones
      w_hat <- as.numeric(w_unnorm / sum(w_unnorm))
      sr_test <- realized_sharpe(w_hat, test_data)
      sr_sum <- sr_sum + sr_test
    }
    scores[h_i] <- sr_sum
  }
  best_idx <- which.max(scores)
  h_cv <- candidate_h[best_idx]
  return(h_cv)
}

