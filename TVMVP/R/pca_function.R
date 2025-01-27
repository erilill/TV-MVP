
#' @export
determine_factors <- function(returns, max_R, g_func = function(N, T) log(N * T) / (N * T)) {
  T <- nrow(returns)
  N <- ncol(returns)

  # Initialize storage
  V <- numeric(max_R)
  penalty <- numeric(max_R)

  # Loop over possible number of factors (R)
  for (R in 1:max_R) {
    residuals <- matrix(NA, nrow = T, ncol = N)
    for (r in 1:T){
      # Step 1: Perform PCA with R factors
      pca_result <- local_pca(returns, r = r, bandwidth = bandwidth, m = R, kernel_func = epanechnikov_kernel)
      X_r <- matrix(0, nrow = T, ncol = N)
      X_r <- sweep(returns, 1, sqrt(pca_result$w_r), `*`)
      Lambda_breve_R <- t((1/T*N)*t(X_r)%*%X_r%*%pca_result$loadings)
      F_breve_R <- solve((Lambda_breve_R)%*%t(Lambda_breve_R))%*%(Lambda_breve_R)%*%returns[r,]

      # Step 2: Compute SSR (Sum of Squared Residuals)
      residuals[r,] <- returns[r,] - t(F_breve_R) %*% (Lambda_breve_R)
      V[R] <- sum(residuals^2) / (N * T)

      penalty[R] <- R * g_func(N, T)
    }
  }
  # Step 4: Determine optimal number of factors
  IC_values <- log(V) + penalty
  optimal_R <- which.min(IC_values)

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
local_pca <- function(returns, r, bandwidth, m, kernel_func){
  T <- nrow(returns)
  p <- ncol(returns)

  k_h <- sapply(1:T, function(t) boundary_kernel(r, t, T, bandwidth, kernel_func))
  X_r <- matrix(0, nrow = T, ncol = p)
  X_r <- sweep(returns, 1, sqrt(k_h), `*`)


  eigen_txr_xr <- eigen((X_r)%*%t(X_r))
  idx <- order(eigen_txr_xr$values, decreasing = TRUE)
  eigvals <- eigen_txr_xr$values[idx]
  eigvecs <- eigen_txr_xr$vectors[, idx]

  F_hat_r <- sqrt(T)*eigvecs[,1:m]

  t_lambda_hat_r <- t(F_hat_r)%*%(X_r)/T

  F_hat <- solve((t_lambda_hat_r)%*%t(t_lambda_hat_r))%*%(t_lambda_hat_r)%*%returns[r,]

  return(list(factors = t(F_hat),
              loadings = t(t_lambda_hat_r),
              w_r = k_h))
}
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

#' Perform Local Principal Component Analysis (PCA) with Optimal Factor Selection
#'
#' This function conducts a comprehensive local PCA on asset returns by first selecting
#' the optimal number of factors using an information criterion and then extracting the
#' corresponding factor scores and loadings for each time period within a specified bandwidth.
#'
#' @param returns A numeric matrix of asset returns with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param bandwidth A numeric value indicating the bandwidth parameter defining the window around each time period.
#' Defaults to the value computed using Silverman's rule of thumb.
#' @param max_factors An integer specifying the maximum number of factors to consider during optimal factor selection. Defaults to \code{10}.
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
  p <- ncol(returns)
  T <- nrow(returns)

  # Initialize storage
  factors <- matrix(NA, nrow = T, ncol = m)
  loadings <- vector("list", T)
  weights_list <- vector("list", T)

  # For each time t, do local PCA
  for (t_i in 1:T) {
    local_result <- local_pca(returns, t_i, bandwidth, m, kernel_func)
    factors[t_i, ] <- local_result$factors
    loadings[[t_i]] <- local_result$loadings
    weights_list[[t_i]] <- local_result$w_r
  }

  return(list(
    factors = factors,    # T x m
    loadings = loadings,  # list of length T, each p x m
    m = m,
    weights = weights_list
  ))
}




