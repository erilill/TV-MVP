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


