sqrt_matrix <- function(A) {
  eig <- eigen(A)
  eig_values_sqrt <- diag(sqrt(abs(eig$values)))  # Ensure non-negative eigenvalues
  return(eig$vectors %*% eig_values_sqrt %*% t(eig$vectors))
}

#'
compute_sigma_0 <- function(res_set, T, p) {
  sigma_0 <- crossprod(res_set) / T

  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      sigma_0[i, j] <- sigma_0[i, j] * ((1 - 0.01)^(j - i))
      sigma_0[j, i] <- sigma_0[i, j]
    }
  }

  return(sigma_0)
}

#' the function will return the size of obj
#' and it is smart in the sense that it will choose the suitable unit
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
