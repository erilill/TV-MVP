#' @export
cv_spcov_lambda <- function(R, k = 5, tau = 1e-4, lambda_grid = c(seq(0.05, 1.5, by = 0.05), seq(10, 100, by = 10), seq(150, 1000, by = 50), seq(1250, 6000, by = 250))) {
  T <- nrow(R)  # Number of time points (samples)
  p <- ncol(R)  # Number of variables (assets)
  
  # Generate random permutation of indices for k-fold cross-validation
  rand_rank <- sample(T)  
  subnum <- floor(T / k)  # Number of samples per fold
  
  all_d <- numeric(length(lambda_grid))  # Store cross-validation scores
  
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    P <- lambda * (matrix(1, p, p) - diag(p))  # Create penalty matrix
    diff_sum <- 0  # Store accumulated difference for CV scoring
    
    for (j in seq_len(k)) {
      # Split data into training and testing sets
      test_idx <- rand_rank[((j - 1) * subnum + 1):(j * subnum)]
      train_data <- R[-test_idx, ]  # Remove test samples
      test_data <- R[test_idx, ]    # Hold-out set
      
      # Compute covariance matrix on training data
      S_train <- cov(train_data)  # Sample covariance from training set
      
      # Solve sparse covariance estimation using spcov
      Sigma_init <- diag(diag(S_train))  # Initialize with diagonal
      result <- tryCatch({
        spcov(Sigma_init, S_train, lambda = P, step.size = 0.001, tol.outer = 1e-6, n.inner.steps = 200, n.outer.steps = 200, thr.inner = 1e-3)
      }, error = function(e) {
        return(NULL)  # Skip if spcov fails
      })
      
      if (is.null(result)) next  # Skip this fold if estimation fails
      
      sp_cov <- result$Sigma  # Extract estimated covariance matrix
      
      # Compute covariance matrix from test data
      S_test <- cov(test_data)
      
      # Log-determinant scoring function
      log_det_spcov <- tryCatch({
        determinant(sp_cov, logarithm = TRUE)$modulus
      }, error = function(e) NaN)
      
      # Compute the evaluation criterion
      if (is.nan(log_det_spcov) || log_det_spcov == -Inf) {
        Dt <- NaN
      } else {
        Dt <- -log_det_spcov - sum(diag(solve(sp_cov, S_test)))
      }
      
      # Accumulate CV score
      diff_sum <- diff_sum + Dt
    }
    
    # Store the averaged cross-validation score
    all_d[i] <- diff_sum / k

  }
  
  # Select the best lambda that maximizes the cross-validation score
  best_lambda <- lambda_grid[which.max(all_d)]
  
  return(list(optimal_lambda = best_lambda, all_d = all_d))
}