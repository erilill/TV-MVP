#' Compute the Square Root of a Matrix
#'
#' Computes the square root of a symmetric matrix via eigen decomposition.
#' Negative eigenvalues are handled by taking the square root of their absolute values.
#'
#' @param Amat A numeric symmetric matrix.
#'
#' @return A matrix that is the square root of \code{Amat}.
#'
#' @keywords internal
sqrt_matrix <- function(Amat) {
  eig <- eigen(Amat)
  eig_values_sqrt <- diag(sqrt(abs(eig$values)))  # Ensure non-negative eigenvalues
  return(eig$vectors %*% eig_values_sqrt %*% t(eig$vectors))
}
#' Compute Sigma_0 p.93 Su and Wang (2017).
#'
#' @param res_set Residuals.
#' @param iT Number of time periods.
#' @param ip Number of assets.
#' 
#' @return Sigma_0 from page 93 in Su and Wang (2017).
#' 
#' @keywords internal
compute_sigma_0 <- function(res_set, iT, ip) {
  sigma_0 <- crossprod(res_set) / iT

  for (i in 1:(ip - 1)) {
    for (j in (i + 1):ip) {
      sigma_0[i, j] <- sigma_0[i, j] * ((1 - 0.01)^(j - i))
      sigma_0[j, i] <- sigma_0[i, j]
    }
  }

  return(sigma_0)
}

#' the function will return the size of obj
#' and it is smart in the sense that it will choose the suitable unit
#' @param obj Object
get_object_size <- function(obj) {
  size_bytes <- object.size(obj)  # Get object size in bytes

  # Determine the best unit
  units <- c("B", "KB", "MB", "GB")
  unit_index <- min(floor(log(size_bytes, 1024)), length(units) - 1)

  # Convert size
  size_converted <- size_bytes / (1024 ^ unit_index)

  # Print the result
  return(sprintf("%.2f %s", size_converted, units[unit_index + 1]))
}
