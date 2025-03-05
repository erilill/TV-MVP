
#' @export
determine_factors <- function(returns, max_R, bandwidth) {
  iT <- nrow(returns)
  ip <- ncol(returns)

  # Initialize storage
  V <- numeric(max_R)
  penalty <- numeric(max_R)
  IC_values <- numeric(max_R)

  # Loop over possible number of factors (R)
  for (R in 1:max_R) {
    residuals <- matrix(NA, nrow = iT, ncol = ip)
    prev_F = NULL
    for (r in 1:iT){
      # Step 1: Perform PCA with R factors
      pca_result <- try(local_pca(returns, r = r, bandwidth = bandwidth, 
                                       m = R, kernel_func = epanechnikov_kernel, 
                                       prev_F))
      if("try-error" %in% class(pca_result))
      {
        next
      }
                             

      X_r <- matrix(0, nrow = iT, ncol = ip)
      X_r <- sweep(returns, 1, sqrt(pca_result$w_r), `*`)
      scaled_loadings <- sqrt(ip) * sweep(pca_result$loadings, 2, sqrt(colSums(pca_result$loadings^2)), "/")
      Lambda_breve_R <- t((1/(iT*ip))*t(X_r)%*%X_r%*%scaled_loadings)
      F_breve_R <- solve((Lambda_breve_R)%*%t(Lambda_breve_R))%*%(Lambda_breve_R)%*%returns[r,]

      # Step 2: Compute SSR (Sum of Squared Residuals)
      residuals[r,] <- returns[r,] - t(F_breve_R) %*% (Lambda_breve_R)

      prev_F <- pca_result$F_hat_r
    }
    V[R] <- sum(residuals^2) / (ip * iT)
    penalty[R] <- R * ((ip+iT*bandwidth)/(ip*iT*bandwidth))*log((ip*iT*bandwidth)/(ip+iT*bandwidth))
    IC_values[R] <- log(V[R]) + penalty[R]
  }
  # Step 4: Determine optimal number of factors
  optimal_R <- which.min(IC_values)
  #message(sprintf("Optimal number of factors is %s.", optimal_R))
  return(list(optimal_R = optimal_R, IC_values = IC_values))
}
#' Perform Local Principal Component Analysis (PCA)
#'
#' This function conducts a local PCA on asset returns within a specified bandwidth around a
#' given time point. It extracts a defined number of principal components (factors) and
#' corresponding loadings for portfolio optimization.
#'
#' @param returns A numeric matrix of asset returns with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param r An integer specifying the current time period for which to perform PCA.
#' @param bandwidth A numeric value indicating the bandwidth parameter defining the window
#' around time period \code{r}.
#' @param m An integer specifying the number of principal components (factors) to extract.
#' @param kernel_func A function representing the kernel used for weighting. Typically, an
#' Epanechnikov kernel or another boundary kernel function.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{factors_full}}{Numeric matrix. Factor scores for the extracted principal components.}
#'   \item{\code{loadings_full}}{Numeric matrix. Loadings corresponding to the extracted principal components.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes boundary kernel weights \eqn{w_r} for the current time period \code{r} using the
#'   specified \code{kernel_func} and \code{bandwidth}.
#'   \item Normalizes the weights so that they sum to 1.
#'   \item Applies the square root of the weights to the returns matrix to obtain weighted returns.
#'   \item Performs PCA on the weighted returns to extract the first \code{m} principal components.
#'   \item Normalizes the factor scores and extracts the corresponding loadings.
#' }
#' If the number of extracted principal components is less than 1, the function returns \code{NULL}.
#'
#' @examples
#' # Load necessary library
#' library(ggplot2)
#'
#' # Simulate data for 50 assets over 200 time periods
#' set.seed(123)
#' T <- 200
#' p <- 50
#' returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
#'
#' # Define an Epanechnikov kernel function (assuming it's defined elsewhere)
#' epanechnikov_kernel <- function(u) {
#'   ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
#' }
#'
#' # Perform local PCA for time period r = 100 with bandwidth = 10 and m = 5 factors
#' local_pca_result <- local_pca(returns, r = 100, bandwidth = 10, m = 5, kernel_func = epanechnikov_kernel)
#' print(local_pca_result$factors_full)
#' print(local_pca_result$loadings_full)
#'
#' @export
local_pca <- function(returns, r, bandwidth, m, kernel_func, prev_F = NULL) {
  iT <- nrow(returns)
  ip <- ncol(returns)

  # Compute Kernel Weights
  k_h <- sapply(1:iT, function(t) boundary_kernel(r, t, iT, bandwidth, kernel_func))
  X_r <- sweep(returns, 1, sqrt(k_h), `*`)  # Weighted returns

  # Compute Eigen Decomposition
  eigen_txr_xr <- eigen((X_r) %*% t(X_r))
  idx <- order(eigen_txr_xr$values, decreasing = TRUE)
  eigvals <- eigen_txr_xr$values[idx]
  eigvecs <- eigen_txr_xr$vectors[, idx]

  # Enforce Orthonormality of Factors (F_r)
  F_hat_r <- sqrt(iT) * eigvecs[, 1:m, drop = FALSE]  # (T x m)

  # Align eigenvector directions if previous factors exist
  if (!is.null(prev_F)) {
    for (j in 1:m) {
      if (cor(prev_F[, j], F_hat_r[, j]) < 0) {
        F_hat_r[, j] <- -F_hat_r[, j]  # Flip sign for consistency
      }
    }
  }

  # Compute Loadings: Lambda_r
  t_lambda_hat_r <- t(F_hat_r) %*% (X_r) / iT  # (m x p)
  loadings <- t(t_lambda_hat_r)
    
  # Second pass to compute F_r_hat
  part1 <- crossprod(loadings)   # (m x m)
  part2 <- crossprod(loadings, returns[r, ]) # (m)
  F_r_hat <- solve(part1, part2)        # (m)
  
  return(list(
    factors = F_hat_r,  # (T x m)
    f_hat = t(F_r_hat), # (1 x m)
    loadings = loadings,  # (p x m)
    w_r = as.matrix(k_h)
  ))
}
#' Perform Local Principal Component Analysis (PCA) with Optimal Factor Selection
#'
#' This function conducts a comprehensive local PCA on asset returns by first selecting
#' the optimal number of factors using an information criterion and then extracting the
#' corresponding factor scores and loadings for each time period within a specified bandwidth.
#'
#' @param returns A numeric matrix of asset returns with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param bandwidth A numeric value indicating the bandwidth parameter defining the window around each time period.
#' Defaults to the value computed using Silverman's rule of thumb.
#' @param m An integer specifying the maximum number of factors to consider during optimal factor selection. Defaults to \code{10}.
#' @param kernel_func A function representing the kernel used for weighting. Typically, an
#' Epanechnikov kernel or another boundary kernel function. Defaults to \code{epanechnikov_kernel}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{factors_list}}{List of numeric matrices. Each element corresponds to the factor scores
#'   for a specific time period beyond the initial window.}
#'   \item{\code{loadings_list}}{List of numeric matrices. Each element corresponds to the factor loadings
#'   for a specific time period beyond the initial window.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the optimal number of factors (\code{m}) by minimizing an information criterion
#'   using the \code{select_optimal_factors} function.
#'   \item Iterates over each time period beyond the initial window defined by \code{bandwidth}.
#'   \item For each time period \eqn{t}, performs a local PCA using the \code{local_pca} function
#'   to extract factor scores and loadings.
#'   \item Aggregates the factor scores and loadings into \code{factors_list} and \code{loadings_list}, respectively.
#' }
#' The function also provides informative messages regarding the number of factors selected.
#'
#' @examples
#' # Load necessary libraries
#' library(ggplot2)
#'
#' # Simulate data for 50 assets over 200 time periods
#' set.seed(123)
#' T <- 200
#' p <- 50
#' returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
#'
#' # Define an Epanechnikov kernel function (assuming it's defined elsewhere)
#' epanechnikov_kernel <- function(u) {
#'   ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
#' }
#'
#' # Perform local PCA with default parameters
#' pca_results <- localPCA(returns)
#'
#' # Access factor scores and loadings for time period 100
#' factors_time_100 <- pca_results$factors_list[[100]]
#' loadings_time_100 <- pca_results$loadings_list[[100]]
#'
#' # Display the results
#' print(factors_time_100)
#' print(loadings_time_100)
#'
#' @export
localPCA <- function(returns,
                     bandwidth,
                     m,
                     kernel_func = epanechnikov_kernel) {
  ip <- ncol(returns)
  iT <- nrow(returns)

  # Initialize storage
  factors <- vector("list", iT)
  loadings <- vector("list", iT)
  weights_list <- vector("list", iT)
  f_hat <- matrix(NA, nrow=iT, ncol=m)

  prev_F <- NULL

  # For each time t, do local PCA
  for (t_i in 1:iT) {
    local_result <- local_pca(returns, t_i, bandwidth, m, kernel_func, prev_F)
    factors[[t_i]] <- local_result$factors
    loadings[[t_i]] <- local_result$loadings
    weights_list[[t_i]] <- local_result$w_r
    f_hat[t_i,] <- local_result$f_hat

    prev_F <- local_result$factors
  }

  return(list(
    factors = factors,    # T x m
    loadings = loadings,  # list of length T, each p x m
    m = m,
    weights = weights_list,
    f_hat=f_hat
  ))
}




