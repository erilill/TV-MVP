## In this file, we define the methods of the R6 class


TVMVP$set("public", "set_data", function(data) {
  # for set_data function, the argument data should not be missing
  # this is why I use missing(data) for checking without setting any default value of it
  # note that this is different from initialize function where data can be missing with default value NULL
  if(missing(data)){
    cli::cli_alert_warning("The data is empty.")
  } else{
    private$data <- data
    private$iT <- nrow(data)
    private$ip <- ncol(data)
    
    #tmp_size <- get_object_size(private$data)
    tmp_size <- prettyunits::pretty_bytes(object.size(private$data))
    cli::cli_alert_info("Tibble data set {tmp_size} with {private$iT} rows and {private$ip} columns successfully assigned.")
  }
  
  invisible(self)  # Enables chaining
})


TVMVP$set("public", "get_data", function() {
  if(is.null(private$data)){
    cli::cli_alert_warning("The data is empty.")
    return(invisible(NULL))
  }
  return(private$data)
  # make a copy of the data and return it
})


p_determine_factors <- function(self, private, max_R, bandwidth) {
  T <- private$iT
  N <- private$iN
  
  # Initialize storage
  V <- numeric(max_R)
  penalty <- numeric(max_R)
  IC_values <- numeric(max_R)
  
  # Loop over possible number of factors (R)
  for (R in 1:max_R) {
    residuals <- matrix(NA, nrow = T, ncol = N)
    prev_F = NULL
    for (r in 1:T){
      # Step 1: Perform PCA with R factors
      pca_result <- try(local_pca(returns, r = r, bandwidth = bandwidth, 
                                  m = R, kernel_func = epanechnikov_kernel, 
                                  prev_F))
      if("try-error" %in% class(pca_result))
      {
        next
      }
      
      
      X_r <- matrix(0, nrow = T, ncol = N)
      X_r <- sweep(returns, 1, sqrt(pca_result$w_r), `*`)
      scaled_loadings <- sqrt(N) * sweep(pca_result$loadings, 2, sqrt(colSums(pca_result$loadings^2)), "/")
      Lambda_breve_R <- t((1/(T*N))*t(X_r)%*%X_r%*%scaled_loadings)
      F_breve_R <- solve((Lambda_breve_R)%*%t(Lambda_breve_R))%*%(Lambda_breve_R)%*%returns[r,]
      
      # Step 2: Compute SSR (Sum of Squared Residuals)
      residuals[r,] <- returns[r,] - t(F_breve_R) %*% (Lambda_breve_R)
      
      prev_F <- pca_result$F_hat_r
    }
    V[R] <- sum(residuals^2) / (N * T)
    penalty[R] <- R * ((N+T*bandwidth)/(N*T*bandwidth))*log((N*T*bandwidth)/(N+T*bandwidth))
    IC_values[R] <- log(V[R]) + penalty[R]
  }
  # Step 4: Determine optimal number of factors
  optimal_R <- which.min(IC_values)
  #message(sprintf("Optimal number of factors is %s.", optimal_R))
  return(list(optimal_R = optimal_R, IC_values = IC_values))
}