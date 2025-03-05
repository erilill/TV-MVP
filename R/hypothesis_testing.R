#' Compute \eqn{M_{\hat{}}} Statistic for Covariance Time-Variation Hypothesis Testing
#'
#' This function calculates the \eqn{M_{\hat{}}} statistic, which measures the average squared
#' discrepancy between local and global factor models across all assets and time periods.
#' It quantifies the difference between locally estimated factors/loadings and their global
#' counterparts.
#'
#' @param local_factors A list where each element is a numeric matrix representing the
#' local factor scores for a specific time period. Each matrix should have \eqn{T} rows
#' (time periods) and \eqn{m} columns (factors).
#' @param global_factors A numeric matrix of global factor scores with \eqn{T} rows
#' (time periods) and \eqn{m} columns (factors).
#' @param local_loadings A list where each element is a numeric matrix representing the
#' local factor loadings for a specific time period. Each matrix should have \eqn{N}
#' rows (assets) and \eqn{m} columns (factors).
#' @param global_loadings A numeric matrix of global factor loadings with \eqn{N} rows
#' (assets) and \eqn{m} columns (factors).
#' @param T An integer specifying the number of time periods.
#' @param N An integer specifying the number of assets.
#' @param m An integer specifying the number of factors.
#'
#' @return A numeric scalar \eqn{M_{\hat{}}} representing the average squared discrepancy
#' between local and global factor models across all assets and time periods.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Initializes the \eqn{M_{\hat{}}} statistic to zero.
#'   \item If the number of factors \eqn{m} is equal to one, it ensures that
#'   \code{global_loadings} and \code{global_factors} are treated as matrices.
#'   \item Iterates over each asset \eqn{i = 1} to \eqn{N} and each time period \eqn{t = 1} to \eqn{T}.
#'   \item For each asset and time period, computes:
#'   \itemize{
#'     \item \code{common_H1}: The dot product of the local loadings and local factors.
#'     \item \code{common_H0}: The dot product of the global loadings and global factors.
#'     \item The squared difference \eqn{(common\_H1 - common\_H0)^2} and adds it to \eqn{M_{\hat{}}}.
#'   }
#'   \item After all iterations, normalizes \eqn{M_{\hat{}}} by dividing by the product of \eqn{N} and \eqn{T}.
#' }
#'
#' @examples
#' # Example parameters
#' T <- 100  # Number of time periods
#' N <- 50   # Number of assets
#' m <- 3    # Number of factors
#'
#' # Simulate local factors and loadings
#' local_factors <- lapply(1:T, function(t) matrix(rnorm(m), nrow=1))
#' local_loadings <- lapply(1:T, function(t) matrix(runif(N * m), nrow=N, ncol=m))
#'
#' # Simulate global factors and loadings
#' global_factors <- matrix(rnorm(T * m), nrow=T, ncol=m)
#' global_loadings <- matrix(runif(N * m), nrow=N, ncol=m)
#'
#' # Compute M_hat
#' M_hat <- compute_M_hat(local_factors, global_factors, local_loadings, global_loadings, T, N, m)
#' print(M_hat)
#'
#' @export
compute_M_hat <- function(local_factors, global_factors, local_loadings, global_loadings, iT, ip, m) {
  M_hat <- 0
  if (m == 1){
    global_loadings <- matrix(global_loadings)
    global_factors <- matrix(global_factors)
    local_factors <- matrix(local_factors)
  }
  for (i in 1:ip) {
    for (t in 1:iT) {
      common_H1 <- (local_loadings[[t]][i,]) %*% local_factors[t,]
      common_H0 <- t(global_loadings[i,]) %*% global_factors[t,]
      M_hat <- M_hat + (common_H1 - common_H0)^2
    }
  }
  M_hat <- M_hat / (ip * iT)
  return(M_hat)
}
#' Compute \eqn{B_{pT}} Statistic for Covariance Time-Variation Hypothesis Testing
#'
#' This function calculates the \eqn{B_{pT}} statistic, which is part of the hypothesis
#' testing procedure to determine whether the covariance matrix of asset returns is time-varying.
#' It incorporates kernel-weighted local and global factor interactions along with residuals.
#'
#' @param local_factors A list where each element is a numeric matrix representing the
#' local factor scores for a specific time period. Each matrix should have \eqn{T} rows
#' (time periods) and \eqn{m} columns (factors).
#' @param global_factors A numeric matrix of global factor scores with \eqn{T} rows
#' (time periods) and \eqn{m} columns (factors).
#' @param residuals A numeric matrix of residuals with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param h A numeric value indicating the bandwidth parameter for the kernel function.
#' @param T An integer specifying the number of time periods.
#' @param p An integer specifying the number of assets.
#' @param kernel_func A function representing the kernel used for weighting. Typically, an
#' Epanechnikov kernel or another boundary kernel function.
#'
#' @return A numeric scalar \eqn{B_{pT}} representing the computed statistic based on
#' kernel-weighted factor interactions and residuals.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the sum of squared residuals for each time period \eqn{s}.
#'   \item Constructs the kernel matrix \eqn{K[s,t]} by applying the \code{boundary_kernel}
#'   function to each pair of time periods \eqn{(s,t)}.
#'   \item Calculates the local dot-product matrix \eqn{L[s,t]} as the dot product between
#'   the local factors at time \eqn{s} and \eqn{t}.
#'   \item Computes the global dot-product matrix \eqn{G[s,t]} as the dot product between
#'   the global factors at time \eqn{s} and \eqn{t}.
#'   \item Computes the element-wise squared difference between \eqn{K * L} and \eqn{G},
#'   multiplies it by the residuals, and sums over all \eqn{s,t}.
#'   \item Scales the aggregated value by \eqn{\frac{\sqrt{h}}{T^2 \sqrt{p}}} to obtain \eqn{B_{pT}}.
#' }
#'
#' @examples
#' # Example parameters
#' T <- 100  # Number of time periods
#' p <- 50   # Number of assets
#' m <- 3    # Number of factors
#' h <- 0.1  # Bandwidth parameter
#'
#' # Simulate local factors and loadings
#' local_factors <- lapply(1:T, function(t) matrix(rnorm(m), nrow=1))
#' local_loadings <- lapply(1:T, function(t) matrix(runif(p * m), nrow=p, ncol=m))
#'
#' # Simulate global factors
#' global_factors <- matrix(rnorm(T * m), nrow=T, ncol=m)
#'
#' # Simulate residuals
#' residuals <- matrix(rnorm(T * p, mean=0, sd=0.01), nrow=T, ncol=p)
#'
#' # Define an Epanechnikov kernel function (assuming it's defined elsewhere)
#' epanechnikov_kernel <- function(u) {
#'   ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
#' }
#'
#' # Compute B_pT
#' B_pT <- compute_B_pT(local_factors, global_factors, residuals, h, T, p, epanechnikov_kernel)
#' print(B_pT)
#'
#' @export
compute_B_pT <- function(local_factors, global_factors, residuals, h, iT, ip, kernel_func) {
  res2 <- rowSums(residuals^2)
  K <- outer(1:iT, 1:iT, Vectorize(function(s, t) boundary_kernel(s, t, iT, h, kernel_func)))
  L <- local_factors %*% t(local_factors)
  G <- global_factors %*% t(global_factors)
  D <- (K * L) - G
  D2 <- D^2
  val <- sum(D2 * res2[row(D2)])
  B_pT <- (sqrt(h) / (iT^2 * sqrt(ip))) * val
  
  return(B_pT)
}
#' Compute \eqn{V_{pT}} Statistic for Covariance Time-Variation Hypothesis Testing
#'
#' This function calculates the \eqn{V_{pT}} statistic, which is part of the hypothesis
#' testing procedure to determine whether the covariance matrix of asset returns is time-varying.
#' It incorporates kernel-weighted factor interactions and residual correlations.
#'
#' @param local_factors A list where each element is a numeric matrix representing the
#' local factor scores for a specific time period. Each matrix should have \eqn{T} rows
#' (time periods) and \eqn{m} columns (factors).
#' @param residuals A numeric matrix of residuals with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param h A numeric value indicating the bandwidth parameter for the kernel function.
#' @param T An integer specifying the number of time periods.
#' @param p An integer specifying the number of assets.
#' @param factor_cov A numeric covariance matrix of the factors with dimensions \eqn{m} x \eqn{m}.
#' @param kernel_func A function representing the kernel used for weighting. Typically, an
#' Epanechnikov kernel or another boundary kernel function.
#'
#' @return A numeric scalar \eqn{V_{pT}} representing the computed statistic based on
#' kernel-weighted factor interactions and residual correlations.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Iterates over each pair of time periods \eqn{(s, r)} where \eqn{s < r}.
#'   \item Computes the two-fold convolution kernel value \eqn{\bar{K}_{sr}} using the
#'   \code{two_fold_convolution_kernel} function.
#'   \item Calculates the squared dot product of local factors weighted by the factor covariance
#'   matrix.
#'   \item Computes the squared dot product of residuals between time periods \eqn{s} and \eqn{r}.
#'   \item Aggregates these values across all relevant time period pairs and scales by
#'   \eqn{\frac{2}{T^2 \times p \times h}} to obtain \eqn{V_{pT}}.
#' }
#'
#' @examples
#' # Example parameters
#' T <- 100  # Number of time periods
#' p <- 50   # Number of assets
#' m <- 3    # Number of factors
#' h <- 0.1  # Bandwidth parameter
#'
#' # Simulate local factors
#' local_factors <- lapply(1:T, function(t) matrix(rnorm(m), nrow=1))
#'
#' # Simulate residuals
#' residuals <- matrix(rnorm(T * p, mean=0, sd=0.01), nrow=T, ncol=p)
#'
#' # Simulate factor covariance matrix
#' factor_cov <- matrix(c(1, 0.5, 0.3,
#'                        0.5, 1, 0.4,
#'                        0.3, 0.4, 1), nrow=3, byrow=TRUE)
#'
#' # Define an Epanechnikov kernel function (assuming it's defined elsewhere)
#' epanechnikov_kernel <- function(u) {
#'   ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
#' }
#'
#' # Compute V_pT
#' V_pT <- compute_V_pT(local_factors, residuals, h, T, p, factor_cov, epanechnikov_kernel)
#' print(V_pT)
#'
#' @export
compute_V_pT <- function(local_factors, residuals, h, iT, ip, kernel_func) {
  V_pT <- 0
  for (s in 1:(iT - 1)) {
    for (r in (s + 1):iT) {
      k_bar_sr <- two_fold_convolution_kernel((s - r) / (iT * h), kernel_func)
      term <- k_bar_sr^2 * (t(local_factors[s, ]) %*% (t(local_factors)%*%local_factors*(1/iT)) %*% local_factors[r, ])^2
      V_pT <- V_pT + term * (t(residuals[s, ]) %*% residuals[r, ])^2
    }
  }
  V_pT <- (2 / (iT^2 * ip * h)) * V_pT
  return(V_pT)
}
#' Compute \eqn{J_{pT}} Statistic for Covariance Time-Variation Hypothesis Testing
#'
#' This function calculates the \eqn{J_{pT}} statistic, which is used to test the null
#' hypothesis that the covariance matrix of asset returns is not time-varying. The statistic
#' compares the scaled \eqn{M_{\hat{}}} and \eqn{B_{pT}} against the square root of \eqn{V_{pT}}.
#'
#' @param B_pT A numeric scalar representing the \eqn{B_{pT}} statistic.
#' @param V_pT A numeric scalar representing the \eqn{V_{pT}} statistic.
#' @param M_hat A numeric scalar representing the \eqn{M_{\hat{}}} statistic.
#' @param T An integer specifying the number of time periods.
#' @param p An integer specifying the number of assets.
#' @param h A numeric value indicating the bandwidth parameter used in kernel functions.
#'
#' @return A numeric scalar \eqn{J_{pT}} representing the computed test statistic.
#'
#' @details
#' The \eqn{J_{pT}} statistic is computed using the formula:
#' \deqn{
#' J_{pT} = \frac{T \sqrt{p} \sqrt{h} M_{\hat{}} - B_{pT}}{\sqrt{V_{pT}}}
#' }
#'
#' This statistic is used to assess whether there is significant evidence to reject the null
#' hypothesis of time-invariant covariance matrices. Typically, if \eqn{J_{pT}} exceeds
#' the critical value (e.g., 1.96 for a 5\% significance level), the null hypothesis is rejected.
#'
#' @examples
#' # Example values
#' B_pT <- 0.5
#' V_pT <- 0.1
#' M_hat <- 0.3
#' T <- 100
#' p <- 50
#' h <- 0.1
#'
#' # Compute J_pT
#' J_pT <- compute_J_pT(B_pT, V_pT, M_hat, T, p, h)
#' print(J_pT)
#'
#' @export
compute_J_pT <- function(B_pT, V_pT, M_hat, iT, ip, h) {
  J_pT <- (iT * sqrt(ip) * sqrt(h) * M_hat - B_pT) / sqrt(V_pT)
  return(J_pT)
}
#' Perform Hypothesis Test for Time-Varying Covariance Matrix
#'
#' This function conducts a hypothesis test to determine whether the covariance matrix of
#' asset returns is time-varying. It computes relevant test statistics and outputs the
#' test result based on the \eqn{J_{pT}} statistic.
#'
#' @param local_factors A list where each element is a numeric matrix representing the
#' local factor scores for a specific time period. Each matrix should have \eqn{m} columns
#' (factors) and as many rows as time periods.
#' @param global_factors A numeric matrix of global factor scores with \eqn{T} rows
#' (time periods) and \eqn{m} columns (factors).
#' @param local_loadings A list where each element is a numeric matrix representing the
#' local factor loadings for a specific time period. Each matrix should have \eqn{p}
#' rows (assets) and \eqn{m} columns (factors).
#' @param global_loadings A numeric matrix of global factor loadings with \eqn{p} rows
#' (assets) and \eqn{m} columns (factors).
#' @param residuals A numeric matrix of residuals with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#' @param kernel_func A function representing the kernel used for weighting. Typically, an
#' Epanechnikov kernel or another boundary kernel function. Defaults to \code{epanechnikov_kernel}.
#'
#' @return A numeric scalar \eqn{J_{pT}} representing the computed test statistic.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the factor covariance matrix by aggregating the local factors across all time periods.
#'   \item Determines the optimal bandwidth parameter \code{h} using Silverman's rule of thumb via the \code{silverman} function.
#'   \item Identifies the number of factors \eqn{m} based on the dimensionality of the local factors.
#'   \item Computes the \code{M_hat} statistic using \code{compute_M_hat}.
#'   \item Computes the \code{B_pT} statistic using \code{compute_B_pT}.
#'   \item Computes the \code{V_pT} statistic using \code{compute_V_pT}.
#'   \item Calculates the \code{J_pT} statistic using \code{compute_J_pT}.
#'   \item Outputs a message indicating whether there is evidence that the covariance is time-varying based on the \eqn{J_{pT}} value.
#' }
#'
#' The hypothesis test follows:
#' \itemize{
#'   \item \strong{Null Hypothesis (\eqn{H_0}):} Covariance matrix is time-invariant.
#'   \item \strong{Alternative Hypothesis (\eqn{H_1}):} Covariance matrix is time-varying.
#' }
#'
#' The test uses the standard normal critical value of 1.96 at the 5\% significance level.
#'
#' @examples
#' # Example parameters
#' T <- 100  # Number of time periods
#' p <- 50   # Number of assets
#' m <- 3    # Number of factors
#'
#' # Simulate local factors and loadings
#' local_factors <- lapply(1:T, function(t) matrix(rnorm(m), nrow=1))
#' local_loadings <- lapply(1:T, function(t) matrix(runif(p * m), nrow=p, ncol=m))
#'
#' # Simulate global factors and loadings
#' global_factors <- matrix(rnorm(T * m), nrow=T, ncol=m)
#' global_loadings <- matrix(runif(p * m), nrow=p, ncol=m)
#'
#' # Simulate residuals
#' residuals <- matrix(rnorm(T * p, mean=0, sd=0.01), nrow=T, ncol=p)
#'
#' # Define an Epanechnikov kernel function (assuming it's defined elsewhere)
#' epanechnikov_kernel <- function(u) {
#'   ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
#' }
#'
#' # Perform hypothesis test
#' J_pT_value <- hyptest1(local_factors, global_factors, local_loadings, global_loadings,
#'                        residuals, kernel_func = epanechnikov_kernel)
#' print(J_pT_value)
#'
#' @export
hyptest1 <- function(returns, m, B = 199, kernel_func = epanechnikov_kernel) {
  # Standardize returns
  returns <- scale(returns)
  
  iT <- nrow(returns)
  ip <- ncol(returns)
  h <- silverman(returns)
  
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
  J_pT_bootstrap <- sapply(1:B, function(b) {
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
  
  return(list(J_NT = J_pT, p_value = p_value, J_pT_bootstrap = J_pT_bootstrap))
}
