#' Compute Sum of Squared Residuals (V_m) for Portfolio Optimization
#'
#' This function calculates the sum of squared residuals (\code{V_m}) for portfolio optimization
#' using a rolling window approach. For each time period, it applies boundary kernel weights,
#' performs Principal Component Analysis (PCA) to extract factors, computes residuals, and
#' aggregates the sum of squared residuals across all time periods.
#'
#' @param returns A numeric matrix of asset returns with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param m An integer specifying the number of principal components (factors) to extract.
#' @param kernel_func A function representing the kernel used for weighting. Typically, an
#' Epanechnikov kernel or another boundary kernel function.
#' @param bandwidth A numeric value indicating the bandwidth parameter for the kernel function.
#'
#' @return A numeric scalar representing the sum of squared residuals (\code{V_m}) normalized by
#' the product of the number of assets and time periods (\eqn{p \times T}).
#'
#' @details
#' The function performs the following steps for each time period \eqn{x = 1} to \eqn{T}:
#' \enumerate{
#'   \item Computes boundary kernel weights \eqn{w_x} using the specified \code{kernel_func} and \code{bandwidth}.
#'   \item Normalizes the weights so that they sum to 1.
#'   \item Applies the square root of the weights to the returns matrix.
#'   \item Performs PCA on the weighted returns to extract the first \code{m} principal components.
#'   \item Computes the fitted returns and residuals by subtracting the modeled returns from the actual returns.
#'   \item Aggregates the sum of squared residuals across all time periods.
#' }
#' Finally, \code{V_m} is computed by dividing the total sum of squared residuals by \eqn{p \times T}.
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
#' # Compute V_m for m = 5, bandwidth = 0.1
#' V_m <- compute_V_m(returns, m = 5, kernel_func = epanechnikov_kernel, bandwidth = 0.1)
#' print(V_m)
#'
#' @export
compute_V_m <- function(returns, m, kernel_func, bandwidth) {
  p <- ncol(returns)
  T <- nrow(returns)
  total_ssr <- 0  # Sum of squared residuals

  for (x in 1:T) {
    # Compute boundary kernel weights for time x
    w_x <- sapply(1:T, function(t) boundary_kernel(x, t, T, bandwidth, kernel_func))
    w_x <- w_x / sum(w_x)  # Normalize weights

    # Apply weights to returns
    sqrt_w_x <- sqrt(w_x)
    weighted_returns <- sweep(returns, 1, sqrt_w_x, `*`)

    # Perform PCA
    pca_result <- prcomp(weighted_returns, center = FALSE, scale. = FALSE)
    num_pcs <- min(m, ncol(pca_result$x))
    if (num_pcs < 1) next  # Skip if no PCs are found

    # Extract factor scores and loadings
    Fhat <- pca_result$x[, 1:num_pcs, drop = FALSE] / sqrt(T)  # Normalize
    loadings_hat <- pca_result$rotation[, 1:num_pcs, drop = FALSE]

    # Compute fitted values and residuals
    fitted <- Fhat %*% t(loadings_hat)
    Resid_x <- returns - fitted
    total_ssr <- total_ssr + sum(Resid_x^2)
  }

  # Compute V_m
  V_m <- total_ssr / (p * T)
  return(V_m)
}
#' Select Optimal Number of Factors Using Information Criterion
#'
#' This function determines the optimal number of principal components (factors) to retain
#' in a portfolio optimization context by minimizing an information criterion. The criterion
#' balances model fit and complexity to prevent overfitting.
#'
#' @param returns A numeric matrix of asset returns with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param max_factors An integer specifying the maximum number of factors to consider.
#' @param T_h A numeric value representing the effective window size, typically derived from the bandwidth parameter.
#' @param kernel_func A function representing the kernel used for weighting. Typically, an
#' Epanechnikov kernel or another boundary kernel function.
#' @param bandwidth A numeric value indicating the bandwidth parameter for the kernel function.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{optimal_m}}{Integer. The optimal number of factors selected based on the information criterion.}
#'   \item{\code{IC_values}}{Numeric vector. Information Criterion values for each number of factors from 1 to \code{max_factors}.}
#'   \item{\code{V_m_values}}{Numeric vector. \code{V_m} values corresponding to each number of factors.}
#'   \item{\code{penalty_values}}{Numeric vector. Penalty terms added to the information criterion for each number of factors.}
#' }
#'
#' @details
#' The function iterates over the number of factors \eqn{m} from 1 to \code{max_factors} and performs the following:
#' \enumerate{
#'   \item Computes \code{V_m} using the \code{compute_V_m} function.
#'   \item Calculates a penalty term to account for model complexity.
#'   \item Computes the Information Criterion (IC) as the sum of the logarithm of \code{V_m} and the penalty.
#' }
#' The optimal number of factors \eqn{m} is the one that minimizes the Information Criterion.
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
#' # Select the optimal number of factors with m up to 10 and bandwidth = 0.1
#' selection_results <- select_optimal_factors(
#'   returns = returns,
#'   max_factors = 10,
#'   T_h = T * 0.1,
#'   kernel_func = epanechnikov_kernel,
#'   bandwidth = 0.1
#' )
#' print(selection_results$optimal_m)
#'
#' @export
select_optimal_factors <- function(returns, max_factors, T_h, kernel_func, bandwidth) {
  p <- ncol(returns)
  T <- nrow(returns)

  IC_values <- numeric(max_factors)
  V_m_values <- numeric(max_factors)
  penalty_values <- numeric(max_factors)

  for (m in 1:max_factors) {
    V_m <- compute_V_m(returns, m, kernel_func, bandwidth)
    V_m_values[m] <- V_m

    # Compute penalty
    penalty <- (p + T_h) / (p * T_h) * log((p * T_h) / (p + T_h)) * m
    penalty_values[m] <- penalty

    # Compute Information Criterion (IC)
    IC_values[m] <- log(V_m) + penalty
  }

  optimal_m <- which.min(IC_values)

  return(list(optimal_m = optimal_m, IC_values = IC_values, V_m_values = V_m_values, penalty_values = penalty_values))
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
local_pca <- function(returns, r, bandwidth, m, kernel_func) {
  T <- nrow(returns)
  p <- ncol(returns)

  # Compute boundary kernel weights for time r
  w_r <- sapply(1:T, function(t) boundary_kernel(r, t, T, bandwidth, kernel_func))
  w_r <- w_r / sum(w_r)  # Normalize weights

  # Apply weights to returns
  sqrt_w_r <- sqrt(w_r)
  weighted_returns <- sweep(returns, 1, sqrt_w_r, `*`)

  # Perform PCA
  pca_result <- prcomp(weighted_returns, center = FALSE, scale. = FALSE)

  # Determine actual number of factors
  num_factors <- min(m, ncol(pca_result$x))
  if (num_factors < 1) return(NULL)  # Return NULL if no factors are found

  # Extract factor scores and loadings
  Fhat_all <- pca_result$x[, 1:num_factors, drop = FALSE] / sqrt(T)  # Normalize
  loadings_all <- pca_result$rotation[, 1:num_factors, drop = FALSE]

  return(list(factors_full = Fhat_all, loadings_full = loadings_all))
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
localPCA <- function(returns, bandwidth = silverman(returns), max_factors = 10, kernel_func = epanechnikov_kernel){
  p <- ncol(returns)
  T <- nrow(returns)

  # Select the optimal number of factors
  m_selection <- select_optimal_factors(
    returns = returns,
    max_factors = max_factors,
    T_h = T * bandwidth,
    kernel_func = kernel_func,
    bandwidth = bandwidth
  )

  m <- m_selection$optimal_m

  # Print the number of factors chosen
  message(sprintf("Optimal number of factors selected: %d", m))

  # Initialize lists to store factors and loadings
  factors_list <- vector("list", T)
  loadings_list <- vector("list", T)

  # Perform Local PCA for each time point
  for (t in 1:T) {
    pca_result <- local_pca(returns, t, bandwidth, m, kernel_func)

    # Handle cases where PCA might return NULL
    if (!is.null(pca_result)) {
      factors_list[[t]] <- pca_result$factors_full
      loadings_list[[t]] <- pca_result$loadings_full
    } else {
      factors_list[[t]] <- NA
      loadings_list[[t]] <- NA
    }
  }

  # Return the factors and loadings lists
  return(list(
    factors_list = factors_list,
    loadings_list = loadings_list
  ))
}
