#' @import quadprog
#' 
solve_minvar_portfolio <- function(Sigma) {
  p <- ncol(Sigma)
  ones <- matrix(1, nrow = p, ncol = 1)
  
  # Solve using quadprog
  result <- solve.QP(
    Dmat = as.matrix(Sigma),
    dvec = rep(0, p),
    Amat = ones,
    bvec = 1,
    meq = 1
  )
  
  return(result$solution)  # Extract optimal weights
}

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


