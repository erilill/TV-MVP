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
#' @param iT An integer specifying the number of time periods.
#' @param ip An integer specifying the number of assets.
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
#
#'
#' @keywords internal
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
#' @param iT An integer specifying the number of time periods.
#' @param ip An integer specifying the number of assets.
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
#'
#' @keywords internal
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
#' @param iT An integer specifying the number of time periods.
#' @param ip An integer specifying the number of assets.
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
#'   \eqn{\frac{2}{T^2 × p × h}} to obtain \eqn{V_{pT}}.
#' }
#'
#'
#' @keywords internal
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
#' Test for Time-Varying Covariance via Local PCA and Bootstrap
#'
#' This function performs a hypothesis test for time-varying covariance in asset returns based on Su and Wang (2017).
#' It first standardizes the input returns and then computes a time-varying covariance estimator
#' using a local principal component analysis (Local PCA) approach. The test statistic \eqn{J_{pT}}
#' is computed and its significance is assessed using a bootstrap procedure. The procedure is available either as a stand-alone
#' function or as a method in the `TVMVP` R6 class.
#'
#' @param returns A numeric matrix of asset returns with dimensions \eqn{T × p} (time periods by assets).
#' @param m Integer. The number of factors to extract in the local PCA. See \code{\link{determine_factors}}.
#' @param B Integer. The number of bootstrap replications to perform. Default is 200
#' @param kernel_func Function. A kernel function for weighting observations in the local PCA. Default is \code{epanechnikov_kernel}.
#'
#' @return A list containing:
#' \item{J_NT}{The test statistic \eqn{J_{pT}} computed on the original data.}
#' \item{p_value}{The bootstrap p-value, indicating the significance of time variation in covariance.}
#' \item{J_pT_bootstrap}{A numeric vector of bootstrap test statistics from each replication.}
#'
#' @details
#' Two usage styles:
#' 
#' \preformatted{
#' # Function interface
#' hyptest1(returns, m=2)
#'
#' # R6 method interface
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#' tv$determine_factors(max_m=5)
#' tv$hyptest()
#' tv
#' tv$get_bootstrap()         # prints bootstrap test statistics
#' }
#'    
#' When using the method form, if `m` are omitted,
#' they default to values stored in the object. Results are cached and
#' retrievable via class methods.
#' 
#' The function follows the steps below:
#' 
#' \enumerate{
#'   \item Standardizes the returns.
#'   \item Computes the optimal bandwidth using the Silverman rule.
#'   \item Performs a local PCA on the standardized returns to extract local factors and loadings.
#'   \item Computes a global factor model via singular value decomposition (SVD) to obtain global factors.
#'   \item Calculates residuals by comparing the local PCA reconstruction to the standardized returns.
#'   \item Computes a test statistic \eqn{J_{pT}} based on a function of the residuals and covariance estimates as:
#' 
#'   \deqn{\hat{J}_{pT} = \frac{T p^{1/2} h^{1/2} \hat{M} - \hat{\mathbb{B}}_{pT}}{\sqrt{\hat{\mathbb{V}}_{pT}}},}
#' 
#'   where:
#' 
#'   \deqn{\hat{M} = \frac{1}{pT} \sum_{i=1}^p \sum_{t=1}^T \left(\hat{\lambda}_{it}' \hat{F}_t - \tilde{\lambda}_{i0}' \tilde{F}_t\right),}
#' 
#'   \deqn{\hat{\mathbb{B}}_{pT} = \frac{h^{1/2}}{T^2 p^{1/2}} \sum_{i=1}^p \sum_{t=1}^T \sum_{s=1}^T \left(k_{h,st} \hat{F}_s' \hat{F}_t - \tilde{F}_s' \tilde{F}_t\right)^2 \hat{e}_{is}^2,}
#' 
#'   and
#' 
#'   \deqn{\hat{\mathbb{V}}_{pT} = \frac{2}{p h T^2} \sum_{1\leq s \neq r \leq T} \bar{k}_{sr}^2 \left(\hat{F}_s' \hat{\Sigma}_F \hat{F}_r \right)^2 \left(\hat{e}_r' \hat{e}_s \right)^2.}
#' 
#'   \item A bootstrap procedure is then used to compute the distribution of \eqn{J_{pT}} and derive a p-value.
#' }
#'
#' The function prints a message indicating the strength of evidence for time-varying covariance based on the p-value.
#' 
#' @section References: 
#' Su, L., & Wang, X. (2017). On time-varying factor models: Estimation and testing. Journal of Econometrics, 198(1), 84–101
#'
#' @examples
#' \dontrun{
#' # Simulate some random returns (e.g., 100 periods, 30 assets)
#' set.seed(123)
#' returns <- matrix(rnorm(100*30, mean = 0, sd = 0.02), nrow = 100, ncol = 30)
#'
#' # Test for time-varying covariance using 3 factors and 200 bootstrap replications
#' test_result <- hyptest1(returns, m = 3, B = 200, kernel_func = epanechnikov_kernel)
#'
#' # Print test statistic and p-value
#' print(test_result$J_NT)
#' print(test_result$p_value)
#' 
#' # Or use R6 method interface
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#' tv$determine_factors(max_m=5)
#' tv$hyptest(B = 200, kernel_func = epanechnikov_kernel)
#' tv
#' tv$get_bootstrap()         # prints bootstrap test statistics
#' }
#'
#' @importFrom stats rnorm
#' @export
hyptest1 <- function(returns, m, B = 200, kernel_func = epanechnikov_kernel) {
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
