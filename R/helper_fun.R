#'
find_smallest_matrix <- function(matrix_list) {
  if (length(matrix_list) == 0) {
    stop("The list is empty.")
  }
  # Extract dimensions of all matrices
  dims <- sapply(matrix_list, function(mat) c(nrow(mat), ncol(mat)))
  # Compute total number of elements (rows * cols)
  total_elements <- apply(dims, 2, prod)
  # Find the index of the smallest matrix
  smallest_index <- which.min(total_elements)
  # Return the smallest matrix
  return(dim(matrix_list[[smallest_index]]))
}

#'
handle_cv_bandwidth <- function(returns, m, candidate_h, kernel_func) {
  h <- try(cv_bandwidth(returns, m, candidate_h, kernel_func), silent = TRUE)

  if (inherits(h, "try-error")) {  # Check if an error occurred
    message("Error detected in bandwidth cross-validation. Attempting to lower m...")

    # Extract numeric values from the error message
    error_message <- as.character(h)
    m_eff_candidates <- as.numeric(unlist(regmatches(error_message, gregexpr("[0-9]+", error_message))))

    # Filter valid m_eff values (remove zeros, negatives, and NAs)
    valid_m_eff <- m_eff_candidates[m_eff_candidates > 0 & !is.na(m_eff_candidates)]

    if (length(valid_m_eff) > 0) {
      new_m <- min(valid_m_eff)  # Pick the smallest valid m_eff
      message(sprintf("Switching to m = %d.", new_m))

      return(cv_bandwidth(returns, new_m, candidate_h, kernel_func)$optimal_h)
    } else {
      stop("Failed to extract a valid m_eff from the error message. Bandwidth selection aborted.")
    }
  }

  return(h$optimal_h)
}

#'
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
