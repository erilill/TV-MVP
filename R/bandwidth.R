#' Work in progress
#' Compute Bandwidth Parameter Using Silverman's Rule of Thumb
#'
#' This function calculates the bandwidth parameter for kernel functions using Silverman's rule of thumb,
#' which is commonly used in kernel density estimation to determine an appropriate bandwidth.
#'
#' @param returns A numeric matrix of asset returns with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param T An optional integer specifying the number of time periods. If provided, it overrides
#' \code{returns}.
#' @param p An optional integer specifying the number of assets. If provided, it overrides
#' \code{returns}.
#'
#' @return A numeric value representing the computed bandwidth parameter based on Silverman's rule.
#'
#' @details
#' Silverman's rule of thumb for bandwidth selection is given by:
#' \deqn{bandwidth = \frac{2.35}{\sqrt{12}} \times T^{-0.2} \times p^{-0.1}}
#' where \eqn{T} is the number of time periods and \eqn{p} is the number of assets.
#'
#' If the number of time periods \code{T} and the number of assets \code{p} are not provided,
#' the function extracts these values from the \code{returns} matrix.
#'
#' @examples
#' # Simulate data for 50 assets over 200 time periods
#' set.seed(123)
#' T <- 200
#' p <- 50
#' returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
#'
#' # Compute bandwidth using Silverman's rule of thumb
#' bw <- silverman(returns)
#' print(bw)
#'
#' # Alternatively, provide T and p directly
#' bw_direct <- silverman(returns = NULL, T = 200, p = 50)
#' print(bw_direct)
#'
#' @export
silverman <- function(returns, T=NULL, p=NULL){
  if (!is.null(returns)){
    p <- ncol(returns)
    T <- nrow(returns)
  }
  bandwidth <- (2.35/sqrt(12)) * T^(-0.2) * p^(-0.1)
  return(bandwidth)
}
#'
#' @export
cv_bandwidth <- function(returns, m, candidate_h, kernel_func) {
  P <- ncol(returns)
  T <- nrow(returns)

  if (P / T + 0.1 > 0.7) {
    p_T <- P / T
    test_loop <- round((1 - P / T - 0.1) / 0.1)
  } else {
    p_T <- 0.6
    test_loop <- 3
  }

  testing_set_list <- vector("list", test_loop)
  training_set_list <- vector("list", test_loop)
  loadings_list <- list()

  for (kk in 1:test_loop) {
    start_test <- floor(T * (p_T + 0.1 * kk)) + 1
    end_test <- ifelse(p_T + 0.1 * (kk + 1) > 1, T, floor(T * (p_T + 0.1 * (kk + 1))))

    testing_set_list[[kk]] <- returns[start_test:end_test, ]
    training_set_list[[kk]] <- returns[1:(start_test - 1), ]
  }

  scores <- numeric(length(candidate_h))

  for (h_i in seq_along(candidate_h)) {
    h_val <- candidate_h[h_i]
    sr_sum <- 0

    for (j in seq_len(test_loop)) {

      train_data <- training_set_list[[j]]
      test_data  <- testing_set_list[[j]]

      local_pca_train <- localPCA(train_data, bandwidth = h_val, m = m, kernel_func = kernel_func)
      factor_cov <- cov(local_pca_train$f_hat)
      
      # Check Effective m
      min_m_eff <- (find_smallest_matrix(local_pca_train$loadings)[2])
      if (min_m_eff < m) {
        stop(sprintf("Effective m (%d) is smaller than m (%d).", min_m_eff, m))
      }

      res <- residuals(local_pca_train$f_hat, local_pca_train$loadings, train_data)
      residual_cov <- tryCatch({
        estimate_residual_cov(res)
      }, error = function(e) {
        message("Warning: Singular residual covariance, using identity matrix.")
        diag(P)  # Use identity matrix if estimation fails
      })

      valid_loadings <- Filter(function(x) !anyNA(x), local_pca_train$loadings)
      avg_loadings <- Reduce(`+`, valid_loadings) / length(valid_loadings)

      Sigma_hat <- avg_loadings %*% factor_cov %*% t(avg_loadings) + residual_cov

      Sigma_hat <- Sigma_hat + diag(1e-6, P)

      inv_cov <- solve(Sigma_hat + diag(1e-6, P))  # Regularization
      ones <- rep(1, P)
      w_hat <- as.numeric(inv_cov %*% ones / sum(inv_cov %*% ones))

      real_returns <- rowSums(test_data * w_hat)
      sr_test <- mean(real_returns) / sd(real_returns)
      sr_sum <- sr_sum + sr_test
    }

    scores[h_i] <- sr_sum
  }

  # Select the best bandwidth
  best_idx <- which.max(scores)
  h_cv <- candidate_h[best_idx]

  return(list(optimal_h = h_cv, scores = scores))
}
