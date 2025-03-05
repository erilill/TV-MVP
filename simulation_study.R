####################
# Simulation study #
####################

load_all() # Switch to library(TVMVP) when finished
library(MASS)


# Will mostly be copying Su and Wang DGP1-6, then two simulated based on real data if I can find data

###############
#### DGP's ####
###############

# Factors
generate_factors <- function(T) {
  F1 <- numeric(T)
  F2 <- numeric(T)
  F1[1] <- rnorm(1, mean = 0, sd = sqrt(1 / (1 - 0.6^2)))
  F2[1] <- rnorm(1, mean = 0, sd = sqrt(1 / (1 - 0.3^2)))
  
  for (t in 2:T) {
    F1[t] <- 0.6 * F1[t - 1] + rnorm(1, mean = 0, sd = sqrt(1 - 0.6^2))
    F2[t] <- 0.3 * F2[t - 1] + rnorm(1, mean = 0, sd = sqrt(1 - 0.3^2))
  }
  
  return(cbind(F1, F2))
}

#DGP 1 IID
generate_DGP1 <- function(p, T, F){
  X <- matrix(NA, nrow=T, ncol=p)
  for (t in 1:T){
    Lambda <- MASS::mvrnorm(p, mu=rep(0,2), Sigma = diag(1,2))
    e <- rnorm(p, mean=0, sd=1)
    X[t,] <- t(Lambda %*% F[t,] + e)
    }
  return(X)
}

# DGP 2 heteroskedasticity
generate_DGP2 <- function(p, T, F){
  X <- matrix(NA, nrow=T, ncol=p)
  for (t in 1:T){
    Lambda <- MASS::mvrnorm(p, mu=rep(0,2), Sigma = diag(1,2))
    sigma_i <- runif(p, 0.5, 1.5)
    e <- rnorm(p, mean=0, sd=1)*sigma_i
    X[t,] <- t(Lambda %*% F[t,] + e)
  }
  return(X)
}

# DGP 3 Cross-sectional dependence
generate_DGP3 <- function(p, T, F){
  X <- matrix(NA, nrow=T, ncol=p)
  for (t in 1:T){
    Lambda <- MASS::mvrnorm(p, mu=rep(0,2), Sigma = diag(1,2))
    Sigma_e <- outer(1:p, 1:p, function(i, j) 0.5^abs(i-j))
    e <- MASS::mvrnorm(1, mu=rep(0,p), Sigma=Sigma_e)
    X[t,] <- t(Lambda %*% F[t,] + e)
  }
  return(X)
}

# DGP 4 Single structural brake
generate_DGP4 <- function(p, T, F, b = 2) {
  X <- matrix(NA, nrow=T, ncol=p)
  for (t in 1:T){
    Lambda <- matrix(rnorm(p * 2, mean = 1, sd = 1), ncol = 2)
    if ((T/2+1) <= t & t <= T){
      Lambda <- Lambda + b
    }
    sigma_i <- runif(p, 0.5, 1.5)
    e <- rnorm(p, mean=0, sd=1)*sigma_i
    X[t,] <- t(Lambda %*% F[t,] + e)
  }
  return(X)
}




# DGP 5 Multiple Structural Breaks
generate_DGP5 <- function(p, T, F, b = 2) {
  X <- matrix(NA, nrow=T, ncol=p)
  for (t in 1:T) {
    Lambda1 <- matrix(rnorm(p , mean = 1, sd = 1), ncol = 1)
    Lambda2 <- matrix(rnorm(p , mean = 0, sd = 1), ncol = 1)
    Lambda <- cbind(Lambda1, Lambda2)
    
    if (0.6 * T < t & t <= 0.8 * T) {
      Lambda[,1] <- Lambda[,1] + 0.5 * b
    } else if (0.2 * T < t & t <= 0.4 * T) {
      Lambda[,1] <- Lambda[,1] - 0.5 * b
    }
    e <- rnorm(p, mean=0, sd=1)
    X[t,] <- t(Lambda %*% F[t,] + e)
  }
  return(X)
}

# DGP 6 Smooth Structural Changes I
generate_DGP6 <- function(p, T, b = 2) {
  X <- matrix(NA, nrow=T, ncol=p)
  
  G <- function(x, a, b) exp(-a * (x - b)^-1)  # Smooth transition function
  
  for (t in 1:T) {
    Lambda <- matrix(rnorm(p * 2, mean = 0, sd = 1), ncol = 2)
    for (i in 1:p) {
      Lambda[i, 2] <- b*G(10 * t / T, 2, 5 * i / p + 2)
    }
    e <- rnorm(p, mean=0, sd=1)
    X[t,] <- t(Lambda %*% F[t,] + e)
  }
  return(X)
}
################################################################################
# run simulation for hypothesis test:
library(parallel)

# Set up
set.seed(1337)
T <- 200
p <- 100
R <- 500
m <- 2
B <- 200

# Generate factors
F_t <- generate_factors(T)

##############
#    DGP1    #
##############

# Set up a cluster using all available cores
cl <- makeCluster(detectCores()-1)

# Export required functions and objects to the cluster's environment.
# Make sure to include all functions (generate_DGP1, hyptest1, etc.) 
# and objects (F_t, T, p, m, B) that will be used inside the worker function.
clusterExport(cl, varlist = c("generate_DGP1", "residuals","sqrt_matrix", "hyptest1", "compute_sigma_0", "silverman", "local_pca", "localPCA", "two_fold_convolution_kernel", "boundary_kernel", "compute_B_pT", "compute_M_hat", "compute_V_pT", "epanechnikov_kernel", "F_t", "T", "p", "m", "B", "R"))
clusterEvalQ(cl, {
  library(TVMVP)
  library(MASS)
})

start.time <- Sys.time()
# Run the simulation in parallel using parLapply.
# Each worker generates synthetic data and runs the hypothesis test.
results_list_dgp1 <- parLapply(cl, 1:R, function(r) {
  # Generate synthetic data for replication r
  X_sim <- generate_DGP1(p, T, F_t)
  
  # Run the hypothesis test on the synthetic data
  test_result <- hyptest1(X_sim, m, B)
  
  # Return the relevant results as a list
  list(
    J_NT = test_result$J_NT,
    p_value = test_result$p_value,
    reject_H0 = test_result$p_value < 0.05
  )
})
# Shut down the cluster
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Combine the list of results into a data frame
test_results_dgp1 <- do.call(rbind, lapply(results_list_dgp1, as.data.frame))

# Compute the rejection rate
rejection_rate <- mean(test_results_dgp1$reject_H0)

# Print summary statistics and a histogram of p-values
print(summary(test_results_dgp1$p_value))
hist(test_results_dgp1$p_value, breaks = 20, col = "lightblue", 
     main = "Histogram of p-values", xlab = "p-value", ylab = "Frequency", border = "black")

# Optionally, print the rejection rate
print(sprintf("Rejection rate: %.4f", rejection_rate))

##############
#    DGP2    #
##############

# Set up a cluster using all available cores
cl <- makeCluster(detectCores()-1)

# Export required functions and objects to the cluster's environment.
# Make sure to include all functions (generate_DGP1, hyptest1, etc.) 
# and objects (F_t, T, p, m, B) that will be used inside the worker function.
clusterExport(cl, varlist = c("generate_DGP2", "residuals","sqrt_matrix", "hyptest1", "compute_sigma_0", "silverman", "local_pca", "localPCA", "two_fold_convolution_kernel", "boundary_kernel", "compute_B_pT", "compute_M_hat", "compute_V_pT", "epanechnikov_kernel", "F_t", "T", "p", "m", "B", "R"))
clusterEvalQ(cl, {
  library(TVMVP)
  library(MASS)
})

start.time <- Sys.time()
# Run the simulation in parallel using parLapply.
# Each worker generates synthetic data and runs the hypothesis test.
results_list_dgp2 <- parLapply(cl, 1:R, function(r) {
  # Generate synthetic data for replication r
  X_sim <- generate_DGP2(p, T, F_t)
  
  # Run the hypothesis test on the synthetic data
  test_result <- hyptest1(X_sim, m, B)
  
  # Return the relevant results as a list
  list(
    J_NT = test_result$J_NT,
    p_value = test_result$p_value,
    reject_H0 = test_result$p_value < 0.05
  )
})
# Shut down the cluster
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Combine the list of results into a data frame
test_results_dgp2 <- do.call(rbind, lapply(results_list_dgp2, as.data.frame))

# Compute the rejection rate
rejection_rate <- mean(test_results_dgp2$reject_H0)

# Print summary statistics and a histogram of p-values
print(summary(test_results_dgp2$p_value))
hist(test_results_dgp2$p_value, breaks = 20, col = "lightblue", 
     main = "Histogram of p-values", xlab = "p-value", ylab = "Frequency", border = "black")

# Optionally, print the rejection rate
print(sprintf("Rejection rate: %.4f", rejection_rate))

##############
#    DGP3    #
##############

# Set up a cluster using all available cores
cl <- makeCluster(detectCores()-1)

# Export required functions and objects to the cluster's environment.
# Make sure to include all functions (generate_DGP1, hyptest1, etc.) 
# and objects (F_t, T, p, m, B) that will be used inside the worker function.
clusterExport(cl, varlist = c("generate_DGP3", "residuals","sqrt_matrix", "hyptest1", "compute_sigma_0", "silverman", "local_pca", "localPCA", "two_fold_convolution_kernel", "boundary_kernel", "compute_B_pT", "compute_M_hat", "compute_V_pT", "epanechnikov_kernel", "F_t", "T", "p", "m", "B", "R"))
clusterEvalQ(cl, {
  library(TVMVP)
  library(MASS)
})

start.time <- Sys.time()
# Run the simulation in parallel using parLapply.
# Each worker generates synthetic data and runs the hypothesis test.
results_list_dgp3 <- parLapply(cl, 1:R, function(r) {
  # Generate synthetic data for replication r
  X_sim <- generate_DGP3(p, T, F_t)
  
  # Run the hypothesis test on the synthetic data
  test_result <- hyptest1(X_sim, m, B)
  
  # Return the relevant results as a list
  list(
    J_NT = test_result$J_NT,
    p_value = test_result$p_value,
    reject_H0 = test_result$p_value < 0.05
  )
})
# Shut down the cluster
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Combine the list of results into a data frame
test_results_dgp3 <- do.call(rbind, lapply(results_list_dgp3, as.data.frame))

# Compute the rejection rate
rejection_rate <- mean(test_results_dgp3$reject_H0)

# Print summary statistics and a histogram of p-values
print(summary(test_results_dgp3$p_value))
hist(test_results_dgp3$p_value, breaks = 20, col = "lightblue", 
     main = "Histogram of p-values", xlab = "p-value", ylab = "Frequency", border = "black")

# Optionally, print the rejection rate
print(sprintf("Rejection rate: %.4f", rejection_rate))

##############
#    DGP4    #
##############

# Set up a cluster using all available cores
cl <- makeCluster(detectCores()-1)

# Export required functions and objects to the cluster's environment.
# Make sure to include all functions (generate_DGP1, hyptest1, etc.) 
# and objects (F_t, T, p, m, B) that will be used inside the worker function.
clusterExport(cl, varlist = c("generate_DGP4", "residuals","sqrt_matrix", "hyptest1", "compute_sigma_0", "silverman", "local_pca", "localPCA", "two_fold_convolution_kernel", "boundary_kernel", "compute_B_pT", "compute_M_hat", "compute_V_pT", "epanechnikov_kernel", "F_t", "T", "p", "m", "B", "R"))
clusterEvalQ(cl, {
  library(TVMVP)
  library(MASS)
})

start.time <- Sys.time()
# Run the simulation in parallel using parLapply.
# Each worker generates synthetic data and runs the hypothesis test.
results_list_dgp4 <- parLapply(cl, 1:R, function(r) {
  # Generate synthetic data for replication r
  X_sim <- generate_DGP4(p, T, F_t)
  
  # Run the hypothesis test on the synthetic data
  test_result <- hyptest1(X_sim, m, B)
  
  # Return the relevant results as a list
  list(
    J_NT = test_result$J_NT,
    p_value = test_result$p_value,
    reject_H0 = test_result$p_value < 0.05
  )
})
# Shut down the cluster
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Combine the list of results into a data frame
test_results_dgp4 <- do.call(rbind, lapply(results_list_dgp4, as.data.frame))

# Compute the rejection rate
rejection_rate <- mean(test_results_dgp4$reject_H0)

# Print summary statistics and a histogram of p-values
print(summary(test_results_dgp4$p_value))
hist(test_results_dgp4$p_value, breaks = 20, col = "lightblue", 
     main = "Histogram of p-values", xlab = "p-value", ylab = "Frequency", border = "black")

# Optionally, print the rejection rate
print(sprintf("Rejection rate: %.4f", rejection_rate))

##############
#    DGP5    #
##############

# Set up a cluster using all available cores
cl <- makeCluster(detectCores()-1)

# Export required functions and objects to the cluster's environment.
# Make sure to include all functions (generate_DGP1, hyptest1, etc.) 
# and objects (F_t, T, p, m, B) that will be used inside the worker function.
clusterExport(cl, varlist = c("generate_DGP5", "residuals","sqrt_matrix", "hyptest1", "compute_sigma_0", "silverman", "local_pca", "localPCA", "two_fold_convolution_kernel", "boundary_kernel", "compute_B_pT", "compute_M_hat", "compute_V_pT", "epanechnikov_kernel", "F_t", "T", "p", "m", "B", "R"))
clusterEvalQ(cl, {
  library(TVMVP)
  library(MASS)
})

start.time <- Sys.time()
# Run the simulation in parallel using parLapply.
# Each worker generates synthetic data and runs the hypothesis test.
results_list_dgp5 <- parLapply(cl, 1:R, function(r) {
  # Generate synthetic data for replication r
  X_sim <- generate_DGP5(p, T, F_t)
  
  # Run the hypothesis test on the synthetic data
  test_result <- hyptest1(X_sim, m, B)
  
  # Return the relevant results as a list
  list(
    J_NT = test_result$J_NT,
    p_value = test_result$p_value,
    reject_H0 = test_result$p_value < 0.05
  )
})
# Shut down the cluster
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Combine the list of results into a data frame
test_results_dgp5 <- do.call(rbind, lapply(results_list_dgp5, as.data.frame))

# Compute the rejection rate
rejection_rate <- mean(test_results_dgp5$reject_H0)

# Print summary statistics and a histogram of p-values
print(summary(test_results_dgp5$p_value))
hist(test_results_dgp5$p_value, breaks = 20, col = "lightblue", 
     main = "Histogram of p-values", xlab = "p-value", ylab = "Frequency", border = "black")

# Optionally, print the rejection rate
print(sprintf("Rejection rate: %.4f", rejection_rate))

##############
#    DGP6    #
##############

# Set up a cluster using all available cores
cl <- makeCluster(detectCores()-1)

# Export required functions and objects to the cluster's environment.
# Make sure to include all functions (generate_DGP1, hyptest1, etc.) 
# and objects (F_t, T, p, m, B) that will be used inside the worker function.
clusterExport(cl, varlist = c("generate_DGP6", "residuals","sqrt_matrix", "hyptest1", "compute_sigma_0", "silverman", "local_pca", "localPCA", "two_fold_convolution_kernel", "boundary_kernel", "compute_B_pT", "compute_M_hat", "compute_V_pT", "epanechnikov_kernel", "F_t", "T", "p", "m", "B", "R"))
clusterEvalQ(cl, {
  library(TVMVP)
  library(MASS)
})

start.time <- Sys.time()
# Run the simulation in parallel using parLapply.
# Each worker generates synthetic data and runs the hypothesis test.
results_list_dgp6 <- parLapply(cl, 1:R, function(r) {
  # Generate synthetic data for replication r
  X_sim <- generate_DGP6(p, T, F_t)
  
  # Run the hypothesis test on the synthetic data
  test_result <- hyptest1(X_sim, m, B)
  
  # Return the relevant results as a list
  list(
    J_NT = test_result$J_NT,
    p_value = test_result$p_value,
    reject_H0 = test_result$p_value < 0.05
  )
})
# Shut down the cluster
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Combine the list of results into a data frame
test_results_dgp6 <- do.call(rbind, lapply(results_list_dgp6, as.data.frame))

# Compute the rejection rate
rejection_rate <- mean(test_results_dgp6$reject_H0)

# Print summary statistics and a histogram of p-values
print(summary(test_results_dgp6$p_value))
hist(test_results_dgp6$p_value, breaks = 20, col = "lightblue", 
     main = "Histogram of p-values", xlab = "p-value", ylab = "Frequency", border = "black")

# Optionally, print the rejection rate
print(sprintf("Rejection rate: %.4f", rejection_rate))


################################################################################
#                             Empirical Data                                   #
################################################################################
# Perhaps I will use the empirical data for simulation as well.
# Compute the local factors and then simulate data based on this.

# Read data from excel, skip NA collumns
omx <- read_excel("C:/Users/erikl_xzy542i/Documents/Master_local/Thesis/Data/omx.xlsx", 
                  col_types = c("date", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "skip", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "skip", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "skip", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "skip", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "skip", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "skip", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "skip", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "skip", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "skip", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "skip", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "skip"))
stibor <- read_excel("C:/Users/erikl_xzy542i/Documents/Master_local/Thesis/Data/stibor.xlsx", 
                     col_types = c("date", "numeric"))

library(xts) # Time series format
stibor <- xts(stibor[, -1], order.by = stibor[[1]])/100 # Convert to decimals
omx <- xts(omx[,-1], order.by = omx[[1]])

# Select 100 random stocks (need to decrease dimension)
random100 <- sample(1:381, 100)
test_sample <- omx[, c(random100)]

# Excess returns
returns <- as.matrix(diff(log(test_sample))[-1,])
risk_free <- as.numeric(((1 + stibor)^(1/252) - 1))[-1] # Annualized, correct?

# Data set includes "röda dagar" which need to be removed
# Find indices of rows where all elements are zero
zero_rows <- which(apply(returns, 1, function(x) all(x == 0)))

# Remove "röda dagar"
returns <- returns[-zero_rows,]
risk_free <- risk_free[-zero_rows]

# Constructing function for rolling window comparison, not for package
# Compare TVMVP with MVP constructed with other methods, i.e. sample cov, poet, EMWA, shrinkage, etc

library(parallel)
library(PortfolioMoments)
library(corpcor)
library(POET)
library(glasso)

mega_rol_pred_parallel <- function(returns,
                                   initial_window,   
                                   rebal_period,     
                                   max_factors,
                                   rf = 0,
                                   num_cores = detectCores() - 1) {  # Detect available cores
  # Dimensions
  T <- nrow(returns)
  p <- ncol(returns)
  
  # Define rebalancing dates **before** starting the cluster
  rebalance_dates <- seq(initial_window + 1, T, by = rebal_period)
  RT <- length(rebalance_dates)
  
  if (!is.null(rf)) {
    rf <- rf[initial_window + 1:T]
  }
  
  m <- determine_factors(returns[1:initial_window,], max_factors, silverman(returns[1:initial_window,]))$optimal_R
  
  # Determine optimal number of factors
  # Start parallel cluster
  cl <- makeCluster(num_cores)
  clusterExport(cl, varlist = c("returns", "rebalance_dates", "m", "rebal_period", "p", "rf", "residuals","sqrt_matrix", "compute_sigma_0", "silverman", "local_pca", "localPCA", "two_fold_convolution_kernel", "boundary_kernel", "epanechnikov_kernel", "estimate_residual_cov_poet_local", "adaptive_poet_rho", "determine_factors"), envir = environment())
  clusterEvalQ(cl, {
    library(PortfolioMoments)
    library(corpcor)
    library(POET)
    library(glasso)
  })
  
  results <- parLapply(cl, seq_len(RT),
                       function(l, rebalance_dates, m, rebal_period, p, rf) {
                         reb_t <- rebalance_dates[l]
                         est_data <- returns[1:(reb_t - 1), , drop = FALSE]
                         
                         # 1) 1/N Equal-Weight Portfolio
                         w_equal <- rep(1 / p, p)
                         
                         # 2) Sample Covariance GMVP
                         Sigma_sample <- cov(est_data)
                         inv_sample <- solve(Sigma_sample)
                         w_sample <- as.numeric(inv_sample %*% rep(1, p))
                         w_sample <- w_sample / sum(w_sample)
                         
                         # 3) Shrinkage Covariance GMVP
                         Sigma_shrink <- corpcor::cov.shrink(est_data)
                         inv_shrink <- solve(Sigma_shrink)
                         w_shrink <- as.numeric(inv_shrink %*% rep(1, p))
                         w_shrink <- w_shrink / sum(w_shrink)
                         
                         # 4) EWMA Covariance GMVP
                         lambda <- 0.94
                         Sigma_emwa <- PortfolioMoments::cov_ewma(est_data, lambda = lambda)
                         inv_emwa <- solve(Sigma_emwa)
                         w_emwa <- as.numeric(inv_emwa %*% rep(1, p))
                         w_emwa <- w_emwa / sum(w_emwa)
                         
                         # 5) POET Covariance GMVP
                         poet_res <- POET(t(est_data), m)
                         Sigma_POET <- poet_res$SigmaY
                         inv_POET <- solve(Sigma_POET)
                         w_POET <- as.numeric(inv_POET %*% rep(1, p))
                         w_POET <- w_POET / sum(w_POET)
                         
                         # 6) Glasso Covariance GMVP
                         S <- cov(est_data)
                         glasso_res <- glasso::glasso(S, rho = 0.01)
                         Sigma_glasso <- glasso_res$w
                         inv_glasso <- solve(Sigma_glasso)
                         w_glasso <- as.numeric(inv_glasso %*% rep(1, p))
                         w_glasso <- w_glasso / sum(w_glasso)
                         
                         # 7) TV-MVP
                         local_res <- localPCA(est_data,  silverman(est_data), m, epanechnikov_kernel)
                         Sigma_hat <- estimate_residual_cov_poet_local(local_res, est_data)$total_cov
                         inv_cov <- MASS::ginv(Sigma_hat)
                         ones <- rep(1, p)
                         w_gmv_unnorm <- inv_cov %*% ones
                         w_tvmvp <- as.numeric(w_gmv_unnorm / sum(w_gmv_unnorm))
                         
                         # Define the holding window for returns
                         hold_end <- min(reb_t + rebal_period - 1, T)
                         ret_window <- returns[reb_t:hold_end, , drop = FALSE]
                         
                         list(
                           daily_ret_equal = ret_window %*% w_equal,
                           daily_ret_sample = ret_window %*% w_sample,
                           daily_ret_shrink = ret_window %*% w_shrink,
                           daily_ret_emwa = ret_window %*% w_emwa,
                           daily_ret_POET = ret_window %*% w_POET,
                           daily_ret_glasso = ret_window %*% w_glasso,
                           daily_ret_tvmvp = ret_window %*% w_tvmvp
                         )
                       },
                       rebalance_dates = rebalance_dates, m = m, rebal_period = rebal_period, p = p, rf = rf
  )
  
  
  stopCluster(cl)  # Stop the parallel cluster
  
  # Extract daily returns for each method
  daily_ret_equal   <- unlist(lapply(results, `[[`, "daily_ret_equal"))
  daily_ret_sample  <- unlist(lapply(results, `[[`, "daily_ret_sample"))
  daily_ret_shrink  <- unlist(lapply(results, `[[`, "daily_ret_shrink"))
  daily_ret_emwa    <- unlist(lapply(results, `[[`, "daily_ret_emwa"))
  daily_ret_POET    <- unlist(lapply(results, `[[`, "daily_ret_POET"))
  daily_ret_glasso  <- unlist(lapply(results, `[[`, "daily_ret_glasso"))
  daily_ret_tvmvp   <- unlist(lapply(results, `[[`, "daily_ret_tvmvp"))
  
  # Compute excess returns
  N <- length(daily_ret_equal)
  
  rf_vec <- if (length(rf) == 1) rep(rf, N) else rf[1:N]
  
  er_equal   <- daily_ret_equal   - rf_vec
  er_sample  <- daily_ret_sample  - rf_vec
  er_shrink  <- daily_ret_shrink  - rf_vec
  er_emwa    <- daily_ret_emwa    - rf_vec
  er_POET    <- daily_ret_POET    - rf_vec
  er_glasso  <- daily_ret_glasso  - rf_vec
  er_tvmvp   <- daily_ret_tvmvp   - rf_vec
  
  # Function to compute performance metrics
  compute_metrics <- function(er) {
    CER <- sum(er)
    mu <- mean(er)
    sd_ <- sd(er)
    sharpe <- mu / sd_
    list(CER = CER, mean_excess = mu, sd = sd_, sharpe = sharpe)
  }
  
  stats_equal  <- compute_metrics(er_equal)
  stats_sample <- compute_metrics(er_sample)
  stats_shrink <- compute_metrics(er_shrink)
  stats_emwa   <- compute_metrics(er_emwa)
  stats_POET   <- compute_metrics(er_POET)
  stats_glasso <- compute_metrics(er_glasso)
  stats_tvmvp  <- compute_metrics(er_tvmvp)
  
  methods_stats <- data.frame(
    method = c("1/N", "SampleCov", "ShrinkCov", "EWMA", "POET", "Glasso", "TV-MVP"),
    cumulative_excess = c(stats_equal$CER,
                          stats_sample$CER,
                          stats_shrink$CER,
                          stats_emwa$CER,
                          stats_POET$CER,
                          stats_glasso$CER,
                          stats_tvmvp$CER),
    mean_excess = c(stats_equal$mean_excess,
                    stats_sample$mean_excess,
                    stats_shrink$mean_excess,
                    stats_emwa$mean_excess,
                    stats_POET$mean_excess,
                    stats_glasso$mean_excess,
                    stats_tvmvp$mean_excess),
    sd = c(stats_equal$sd,
           stats_sample$sd,
           stats_shrink$sd,
           stats_emwa$sd,
           stats_POET$sd,
           stats_glasso$sd,
           stats_tvmvp$sd),
    sharpe = c(stats_equal$sharpe,
               stats_sample$sharpe,
               stats_shrink$sharpe,
               stats_emwa$sharpe,
               stats_POET$sharpe,
               stats_glasso$sharpe,
               stats_tvmvp$sharpe)
  )
  
  list(
    daily_returns = list(equal = daily_ret_equal,
                         sample_cov = daily_ret_sample,
                         shrink_cov = daily_ret_shrink,
                         EWMA = daily_ret_emwa,
                         POET = daily_ret_POET,
                         glasso = daily_ret_glasso,
                         tvmvp = daily_ret_tvmvp),
    stats = methods_stats
  )
}

rolling_window_results_day_p <- mega_rol_pred_parallel(returns, 250, 1, rf=risk_free, max_factors = 10)
rolling_window_results_week_p <- mega_rol_pred_parallel(returns, 250, 5, rf=risk_free, max_factors = 10)
rolling_window_results_month_p <- mega_rol_pred_parallel(returns, 250, 21, rf=risk_free, max_factors = 10)






test_rolpred <- rolling_time_varying_mvp(returns, 250, 5, 10, rf=risk_free[251:502])
test_pred <- predict_portfolio(returns, rf=risk_free)












