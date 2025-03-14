#' Determine the Optimal Number of Factors via an Information Criterion
#'
#' This function selects the optimal number of factors for a local principal component
#' analysis (PCA) model of asset returns. It computes an BIC-type information criterion (IC) for each candidate
#' number of factors, based on the sum of squared residuals (SSR) from the PCA reconstruction and a
#' penalty term that increases with the number of factors. The optimal number of factors is chosen as the
#' one that minimizes the IC.
#'
#' @param returns A numeric matrix of asset returns with dimensions \eqn{T \times p}, where \eqn{T} is the number of observations and \eqn{p} is the number of assets.
#' @param max_m Integer. The maximum number of factors to consider.
#' @param bandwidth Numeric. The bandwidth used in the kernel weighting for the local PCA.
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{optimal_R}: The optimal number of factors (an integer) that minimizes the information criterion.
#'   \item \code{IC_values}: A numeric vector of length \code{max_m} containing the information criterion values
#'         for each candidate number of factors.
#' }
#'
#' @details
#' For each candidate number of factors \eqn{R} (from 1 to \code{max_m}), the function:
#'
#' \enumerate{
#'   \item Performs a local PCA on the returns at each time point \eqn{r = 1,\dots,T} using \eqn{R} factors.
#'   \item Computes a reconstruction of the returns and the corresponding residuals:
#'         \deqn{\text{Residual}_r = R_r - F_r \Lambda_r,}
#'         where \eqn{R_r} is the return at time \eqn{r}, and \eqn{F_r} and \(\Lambda_r\) are the local factors and loadings, respectively.
#'   \item Computes the average sum of squared residuals (SSR) as:
#'         \deqn{V(R) = \frac{1}{pT} \sum_{r=1}^{T} \| \text{Residual}_r \|^2.}
#'   \item Adds a penalty term that increases with \eqn{R}:
#'         \deqn{\text{Penalty}(R) = R \times \frac{(p + T \times \text{bandwidth})}{(pT \times \text{bandwidth})} \log\left(\frac{pT \times \text{bandwidth}}{(p + T \times \text{bandwidth})}\right).}
#'   \item The information criterion is defined as:
#'         \deqn{\text{IC}(R) = \log\big(V(R)\big) + \text{Penalty}(R).}
#' }
#'
#' The optimal number of factors is then chosen as the value of \eqn{R} that minimizes \(\text{IC}(R)\).
#'
#' @examples
#' \dontrun{
#' # Simulate a returns matrix with 200 observations and 10 assets
#' set.seed(123)
#' returns <- matrix(rnorm(200 * 10), nrow = 200, ncol = 10)
#'
#' # Determine the optimal number of factors (up to 5) using a specified bandwidth.
#' result <- determine_factors(returns, max_m = 5, bandwidth = 0.2)
#' print(result$optimal_R)
#' print(result$IC_values)
#' }
#'
#' @export
determine_factors <- function(returns, max_m, bandwidth) {
  iT <- nrow(returns)
  ip <- ncol(returns)

  # Initialize storage
  V <- numeric(max_m)
  penalty <- numeric(max_m)
  IC_values <- numeric(max_m)

  # Loop over possible number of factors (R)
  for (mi in 1:max_m) {
    residuals <- matrix(NA, nrow = iT, ncol = ip)
    prev_F = NULL
    for (r in 1:iT){
      # Step 1: Perform PCA with R factors
      pca_result <- try(local_pca(returns, r = r, bandwidth = bandwidth, 
                                       m = mi, kernel_func = epanechnikov_kernel, 
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
    V[mi] <- sum(residuals^2) / (ip * iT)
    penalty[mi] <- mi * ((ip+iT*bandwidth)/(ip*iT*bandwidth))*log((ip*iT*bandwidth)/(ip+iT*bandwidth))
    IC_values[mi] <- log(V[mi]) + penalty[mi]
  }
  # Step 4: Determine optimal number of factors
  optimal_R <- which.min(IC_values)
  #message(sprintf("Optimal number of factors is %s.", optimal_R))
  return(list(optimal_m = optimal_m, IC_values = IC_values))
}
#' Perform Local Principal Component Analysis
#'
#' This function performs a local principal component analysis (PCA) on asset returns,
#' weighted by a specified kernel function. It extracts local factors and loadings from the weighted
#' returns and computes a factor estimate. Optionally, previously estimated factors can
#' be provided to align the new factors' directions.
#'
#' @param returns A numeric matrix of asset returns with dimensions \eqn{T \times p}, where \eqn{T} is the number of time periods and \eqn{p} is the number of assets.
#' @param r Integer. The current time index at which to perform the local PCA.
#' @param bandwidth Numeric. The bandwidth used in the kernel weighting.
#' @param m Integer. The number of factors to extract.
#' @param kernel_func Function. The kernel function used for weighting observations (e.g., \code{epanechnikov_kernel}).
#' @param prev_F Optional. A numeric matrix of previously estimated factors (with dimensions \eqn{T \times m}) used for aligning eigenvector directions. Default is \code{NULL}.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{factors}: A \eqn{T \times m} matrix of local factors estimated from the weighted returns.
#'   \item \code{f_hat}: A \eqn{1 \times m} vector containing the factor estimate for time \eqn{r}.
#'   \item \code{loadings}: A \eqn{p \times m} matrix of factor loadings.
#'   \item \code{w_r}: A numeric vector of kernel weights used in the computation.
#' }
#'
#' @details
#' The function operates in the following steps:
#'
#' \enumerate{
#'   \item **Kernel Weight Computation:**  
#'         For each time point \eqn{t = 1, \dots, T}, the kernel weight is computed using 
#'         \code{boundary_kernel(r, t, T, bandwidth, kernel_func)}. The weighted returns are given by
#'         \deqn{X_r = \text{returns} \odot \sqrt{k_h},}
#'         where \(\odot\) denotes element-wise multiplication and \(k_h\) is the vector of kernel weights.
#'
#'   \item **Eigen Decomposition:**  
#'         The function computes the eigen decomposition of the matrix \eqn{X_r X_r^\top} and orders the eigenvalues in
#'         descending order. The top \eqn{m} eigenvectors are scaled by \(\sqrt{T}\) to form the local factors:
#'         \deqn{\hat{F}_r = \sqrt{T} \, \text{eigvecs}_{1:m}.}
#'
#'   \item **Direction Alignment:**  
#'         If previous factors (\code{prev_F}) are provided, the function aligns the signs of the new factors with the previous ones 
#'         by checking the correlation and flipping the sign if the correlation is negative.
#'
#'   \item **Loadings Computation:**  
#'         The loadings are computed by projecting the weighted returns onto the factors:
#'         \deqn{\Lambda_r = \frac{1}{T} X_r^\top \hat{F}_r,}
#'         where the result is transposed to yield a \eqn{p \times m} matrix.
#'
#'   \item **One-Step-Ahead Factor Estimation:**  
#'         A second pass computes the factor estimate for the current time index \eqn{r} by solving
#'         \deqn{\hat{F}_r = \left(\Lambda_r^\top \Lambda_r\right)^{-1} \Lambda_r^\top R_r,}
#'         where \(R_r\) is the return vector at time \eqn{r}.
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate asset returns for 250 periods and 10 assets:
#' set.seed(123)
#' returns <- matrix(rnorm(250 * 10), nrow = 250, ncol = 10)
#'
#' # Perform local PCA at time index 100 using 4 factors and the Epanechnikov kernel:
#' result <- local_pca(returns, r = 100, bandwidth = 0.2, m = 4, kernel_func = epanechnikov_kernel)
#' str(result)
#' }
#'
#'@keywords internal
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
#' Perform Local PCA Over Time
#'
#' This function performs a local principal component analysis (PCA) on asset returns for each time 
#' period, aggregating the results over time. It calls an internal function \code{local_pca()} at each 
#' time index to extract local factors, loadings, and one-step-ahead factor estimates, and stores these 
#' results in lists. It uses previously computed factors to align the sign of the new factors.
#'
#' @param returns A numeric matrix of asset returns with dimensions \eqn{T \times p}, where \eqn{T} is the 
#' number of time periods and \eqn{p} is the number of assets.
#' @param bandwidth Numeric. The bandwidth parameter used in the kernel weighting for the local PCA.
#' @param m Integer. The number of factors to extract.
#' @param kernel_func Function. The kernel function used for weighting observations. Default is 
#' \code{epanechnikov_kernel}.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{factors}: A list of length \eqn{T}, where each element is a \eqn{T \times m} matrix of local factors.
#'   \item \code{loadings}: A list of length \eqn{T}, where each element is a \eqn{p \times m} matrix of factor loadings.
#'   \item \code{m}: The number of factors extracted.
#'   \item \code{weights}: A list of length \eqn{T}, where each element is a vector of kernel weights used at that time point.
#'   \item \code{f_hat}: A \eqn{T \times m} matrix of one-step-ahead factor estimates.
#' }
#'
#' @details
#' The function processes the input returns over \eqn{T} time periods by iteratively calling the 
#' \code{local_pca()} function. For each time \eqn{t_i}:
#'
#' \enumerate{
#'   \item Kernel weights are computed using the specified \code{kernel_func} and \code{bandwidth}.
#'   \item The returns are weighted by the square root of these kernel weights.
#'   \item An eigen decomposition is performed on the weighted returns' covariance matrix to extract the 
#'         top \eqn{m} eigenvectors, which are scaled by \(\sqrt{T}\) to form the local factors.
#'   \item The signs of the new factors are aligned with those of the previous factors.
#'   \item The factor loadings are computed by projecting the weighted returns onto the local factors, 
#'         normalized by \eqn{T}.
#'   \item A second pass computes a one-step-ahead factor estimate for the current time period.
#' }
#'
#' @examples
#' \dontrun{
#' # Generate simulated returns for 250 periods and 10 assets
#' set.seed(123)
#' returns <- matrix(rnorm(250 * 10), nrow = 250, ncol = 10)
#'
#' # Perform local PCA with a bandwidth of 0.2 and 4 factors using the default Epanechnikov kernel
#' results <- localPCA(returns, bandwidth = 0.2, m = 4)
#' str(results$f_hat)
#' }
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




