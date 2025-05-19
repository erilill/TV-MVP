#' Estimate Local Covariance
#'
#' This internal function computes a time-varying covariance matrix estimate for a given
#' window of asset returns by combining factor-based and sparse residual covariance estimation.
#' It uses results from a local PCA to form residuals and then applies an adaptive thresholding
#' procedure (via \code{adaptive_poet_rho()}) to shrink the residual covariance.
#'
#' @param localPCA_results A list containing the results from local PCA, with components:
#'   \itemize{
#'     \item \code{loadings}: a list where each element is a \eqn{p × m} matrix of factor loadings.
#'     \item \code{f_hat}: a \eqn{T × m} matrix of estimated factors.
#'     \item \code{weights}: a list of kernel weight vectors.
#'   }
#' @param returns A numeric matrix of asset returns with dimensions \eqn{T × p}.
#' @param M0 Integer. The number of observations to leave out between the two sub-samples in the adaptive thresholding procedure. Default is 10.
#' @param rho_grid A numeric vector of candidate shrinkage parameters \eqn{\rho} used in \code{adaptive_poet_rho()}. Default is \code{seq(0.005, 2, length.out = 30)}.
#' @param floor_value A small positive number specifying the lower bound for eigenvalues in the final positive semidefinite repair. Default is \code{1e-12}.
#' @param epsilon2 A small positive tuning parameter for the adaptive thresholding. Default is \code{1e-6}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{best_rho}: The selected shrinkage parameter \eqn{\hat{\rho}_t} for the local residual covariance.
#'     \item \code{residual_cov}: The shrunk residual covariance matrix \eqn{\hat{\Sigma}_e(T)}.
#'     \item \code{total_cov}: The final estimated time-varying covariance matrix \eqn{\Sigma_R(t)}.
#'     \item \code{loadings}: The local factor loadings \eqn{\Lambda_t} from the local PCA.
#'     \item \code{naive_resid_cov}: The raw (unshrunk) residual covariance matrix.
#'   }
#'
#' @details
#' The function follows these steps:
#'
#' \enumerate{
#'   \item **Local Residuals:**
#'         Extract the local loadings \eqn{\Lambda_t} from the last element of \code{localPCA_results\$loadings} and
#'         factors \eqn{\hat{F}} from \code{localPCA_results\$f_hat}. Let \eqn{w_t} denote the corresponding kernel weights.
#'         The local residuals are computed as:
#'         \deqn{U_{\text{local}} = R - F \Lambda_t,}
#'         where \eqn{R} is the returns matrix.
#'
#'   \item **Adaptive Thresholding:**
#'         The function calls \code{adaptive_poet_rho()} on \eqn{U_{\text{local}} }to select an optimal shrinkage parameter
#'         \eqn{\hat{\rho}_t}.
#'
#'   \item **Residual Covariance Estimation:**
#'         The raw residual covariance is computed as:
#'         \deqn{S_{u,\text{raw}} = \frac{1}{T} U_{\text{local}}^\top U_{\text{local}},}
#'         and a threshold is set as:
#'         \deqn{\text{threshold} = \hat{\rho}_t × \text{mean}(|S_{u,\text{raw}}|),}
#'         where the mean is taken over the off-diagonal elements.
#'         Soft-thresholding is then applied to obtain the shrunk residual covariance matrix \eqn{\hat{S}_u}.
#'
#'   \item **Total Covariance Estimation:**
#'         The final covariance matrix is constructed by combining the factor component with the shrunk residual covariance:
#'         \deqn{\Sigma_R(t) = \Lambda_t \left(\frac{F^\top F}{T}\right) \Lambda_t^\top + \hat{S}_u.}
#'
#'   \item **PSD Repair:**
#'         A final positive semidefinite repair is performed by flooring eigenvalues at \code{floor_value} and symmetrizing the matrix.
#' }
#'
#' @examples
#' \dontrun{
#' # Assume localPCA_results is computed via localPCA() and returns is a T x p matrix.
#' cov_est <- estimate_residual_cov_poet_local(localPCA_results, returns)
#' str(cov_est)
#' }
#'
#' @keywords internal
estimate_residual_cov_poet_local <- function(localPCA_results,
                                             returns,
                                             M0 = 10,
                                             rho_grid = seq(0.005, 2, length.out = 30),
                                             floor_value = 1e-12,
                                             epsilon2 = 1e-6) {

  # This function:
  #   1. Form local residuals e_t = R - F * Lambda_t
  #   2. Call adaptive_poet_rho() on those residuals to pick the single best rho_t
  #   3. Compute raw residual covariance, shrink once using rho_t
  #   4. Combine with factor part to get Sigma_R(t)
  #
  # The function returns a list with the final local residual and total covariances.



    # 1. Extract local loadings, factors, and row indices
    Lambda_t <- localPCA_results$loadings[[nrow(returns)]]  # p x m
    F_t <- localPCA_results$f_hat   # T x m
    w_t <- localPCA_results$weights[[nrow(returns)]]

    # 2. residuals
    U_local <- residuals(F_t, localPCA_results$loadings, returns)  # T x p

    # 3. Pick best rho for these local residuals using Chen–Leng–style grouping
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
      residual_cov    = S_u_shrunk,  # Sigma e(T)
      total_cov       = Sigma_R_t,   # Sigma R(T)
      loadings        = Lambda_t,
      naive_resid_cov = S_u_raw
    )


  return(output_list)
}
#' Adaptive Selection of the Shrinkage Parameter \eqn{\rho} for POET
#'
#' This function selects an optimal shrinkage parameter \eqn{\rho} for the residual covariance
#' estimation procedure. It does so by dividing the data into groups and comparing a shrunk covariance
#' matrix (computed on one subsample) to a benchmark covariance (computed on another subsample) using
#' the Frobenius norm. The candidate \eqn{\rho} that minimizes the total squared Frobenius norm difference
#' is selected.
#'
#' @param R A numeric matrix of data (e.g., residuals) with dimensions \eqn{T × p}, where \eqn{T}
#' is the number of observations and \eqn{p} is the number of variables.
#' @param M0 Integer. The number of observations to leave out between two subsamples when forming groups.
#' Default is 10.
#' @param rho_grid A numeric vector of candidate shrinkage parameters \eqn{\rho}. Default is
#' \code{seq(0.001, 2, length.out = 20)}.
#' @param epsilon2 A small positive tuning parameter used as an adjustment in the selection of \eqn{\rho}.
#' Default is \code{1e-6}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{best_rho}: The selected optimal shrinkage parameter \eqn{\hat{\rho}} that minimizes the total
#'   squared Frobenius norm difference.
#'   \item \code{rho_1}: The lower bound for \eqn{\rh} derived from the minimum eigenvalue criteria (adjusted by \code{epsilon2}).
#'   \item \code{min_Fnorm}: The minimum total squared Frobenius norm difference achieved.
#' }
#'
#' @details
#' The function proceeds as follows:
#'
#' \enumerate{
#'   \item The total number of observations \eqn{T} is halved (floored) to define \eqn{T_1} and \eqn{T_2}:
#'         \deqn{T_1 = \left\lfloor \frac{T/2 × \left(1 - 1/\log(T)\right)} \right\rfloor,\quad T_2 = \lfloor T/2 \rfloor - T_1.}
#'
#'   \item The sample is divided into \eqn{\lfloor T/(2M_0) \rfloor} groups (with \eqn{M_0} observations left out in between).
#'
#'   \item For each group, two subsamples are defined:
#'     \itemize{
#'       \item Subsample 1: the first \eqn{T_1} observations of the group.
#'       \item Subsample 2: the last \eqn{T_2} observations of the group, after skipping \eqn{M_0} observations following subsample 1.
#'     }
#'
#'   \item For each group and for a given candidate \eqn{\rho} in \code{rho_grid}, the covariance matrix \eqn{S_1}
#'         is computed from subsample 1, and then shrunk by applying a soft-thresholding function:
#'
#'         \deqn{S_{1,\text{shrunk}} = \text{soft\_threshold}\left(S_1, \rho × \text{mean}\left(|S_1|_{\text{off-diag}}\right)\right).}
#'
#'   \item The function computes the total squared Frobenius norm difference between \eqn{S_{1,\text{shrunk}}}
#'         and the covariance matrix \eqn{S_2} (computed from subsample 2) over all groups.
#'
#'   \item Finally, the function scans across the \code{rho_grid} to select the \eqn{\rho} that minimizes this
#'         total error. In addition, \eqn{\rho_1} is computed as \eqn{\epsilon_2} plus the smallest candidate \eqn{\rho}
#'         for which the smallest eigenvalue of the shrunk covariance is positive.
#' }
#'
#' @examples
#' \dontrun{
#' # Generate a random matrix R (e.g., 100 observations and 10 variables)
#' set.seed(123)
#' R <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#'
#' # Select the optimal shrinkage parameter using the adaptive POET procedure
#' result <- adaptive_poet_rho(R, M0 = 10, rho_grid = seq(0.001, 2, length.out = 20), epsilon2 = 1e-6)
#' print(result)
#' }
#'
#' @importFrom stats cov
#' @keywords internal
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
  return(list(best_rho = best_rho, rho_1 = rho_1, min_Fnorm = min_val))
}


#' Estimate Time-Varying Covariance Matrix Using Local PCA
#'
#' This function estimates a time-varying covariance matrix using local principal component
#' analysis and the soft thresholding for residual shrinkage. By default, only the total
#' covariance matrix is returned. Optionally, the user can retrieve all intermediate
#' components of the estimation process. The procedure is available either as a
#' stand-alone function or as a method in the `TVMVP` R6 class.
#'
#' @param returns A numeric matrix of asset returns with dimensions \eqn{T × p}.
#' @param m The number of factors to use in local PCA.
#' @param bandwidth Optional bandwidth for the local PCA. If not provided, Silverman's rule is used.
#' @param kernel_func The kernel function to use (default is \code{\link{epanechnikov_kernel}}).
#' @param M0 Integer. The number of observations to leave out between the two sub-samples in the adaptive thresholding procedure. Default is 10.
#' @param rho_grid A numeric vector of candidate shrinkage parameters \eqn{\rho} used in \code{\link{adaptive_poet_rho}}. Default is \code{seq(0.005, 2, length.out = 30)}.
#' @param floor_value A small positive number specifying the lower bound for eigenvalues in the final positive semidefinite repair. Default is \code{1e-12}.
#' @param epsilon2 A small positive tuning parameter for the adaptive thresholding. Default is \code{1e-6}.
#' @param full_output Logical; if \code{TRUE}, returns all components of the estimation.
#'
#' @return By default, a covariance matrix. If \code{full_output = TRUE}, a list containing:
#' \itemize{
#'   \item \code{total_cov} – the estimated covariance matrix,
#'   \item \code{residual_cov} – the residual (idiosyncratic) covariance,
#'   \item \code{loadings} – estimated factor loadings,
#'   \item \code{best_rho} – optimal shrinkage parameter,
#'   \item \code{naive_resid_cov} – residual covariance before shrinkage
#' }
#'
#' @details
#' The function estimates a time-varying covariance matrix using Local PCA:
#' \deqn{\hat{\Sigma}_{r,t}=\hat{\Lambda}_t \hat{\Sigma}_F \hat{\Lambda}_t' + \tilde{\Sigma}_e}
#' Where \eqn{\hat{\Lambda}_t} is the factor loadings at time t, \eqn{\hat{\Sigma}_F} is the factor covariance matrix, and \eqn{\tilde{\Sigma}_e} is regularized covariance matrix of the idiosyncratic errors.
#'
#' Two usage styles:
#'
#' \preformatted{
#' # Function interface
#' cov <- time_varying_cov(returns, m = 5)
#'
#' # R6 method interface
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#' tv$determine_factors(max_m = 5)
#' cov <- tv$time_varying_cov()
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' returns <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
#'
#' time_varying_cov(returns = returns, m = 3)
#'
#' # Or using R6 method
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#' tv$time_varying_cov(m=3)
#' }
#'
#' @export
time_varying_cov <- function(returns,
                             m,
                             bandwidth = silverman(returns),
                             kernel_func = epanechnikov_kernel,
                             M0 = 10,
                             rho_grid = seq(0.005, 2, length.out = 30),
                             floor_value = 1e-12,
                             epsilon2 = 1e-6,
                             full_output = FALSE) {
  # Step 1: Local PCA
  local_res <- localPCA(
    returns     = returns,
    bandwidth   = bandwidth,
    m           = m,
    kernel_func = kernel_func
  )

  # Step 2: Residual covariance estimation with POET
  res <- estimate_residual_cov_poet_local(
    localPCA_results = local_res,
    returns           = returns,
    M0                = M0,
    rho_grid          = rho_grid,
    floor_value       = floor_value,
    epsilon2          = epsilon2
  )

  if (full_output) {
    return(res)
  } else {
    return(res$total_cov)
  }
}
