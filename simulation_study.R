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
generate_DGP4 <- function(p, T, F, b = 1) {
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
generate_DGP5 <- function(p, T, F, b = 1) {
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
generate_DGP6 <- function(p, T, b = 1) {
  X <- matrix(NA, nrow=T, ncol=p)
  
  G <- function(x, a, b) exp(-a * (x - b)^2)  # Smooth transition function
  
  for (t in 1:T) {
    Lambda <- matrix(rnorm(p * 2, mean = 0, sd = 1), ncol = 2)
    for (i in 1:p) {
      Lambda[i, 2] <- b * G(10 * t / T, 2, 5 * i / p + 2)
    }
    e <- rnorm(p, mean=0, sd=1)
    X[t,] <- t(Lambda %*% F[t,] + e)
  }
  return(X)
}
################################################################################
# run simulation for hypothesis test:


# Set up
set.seed(1337)
T <- 200
p <- 100
R <- 500
m <- 2
B <- 199

# Storage for results
test_results_dgp1 <- data.frame(J_NT = numeric(R), p_value = numeric(R), reject_H0 = logical(R))

# Generate factors
F_t <- generate_factors(T)

for (r in 1:R) {
  # Generate synthetic data using DGP1
  X_sim <- generate_DGP1(p, T, F_t)
  
  # Run the hypothesis test
  test_result <- hyptest1(X_sim, m, B)
  
  # Store results
  test_results_dgp1$J_NT[r] <- test_result$J_NT
  test_results_dgp1$p_value[r] <- test_result$p_value
  test_results_dgp1$reject_H0[r] <- test_result$p_value < 0.05  # Reject H0 if p-value < 0.05
}

# Compute rejection rate
rejection_rate <- mean(test_results$reject_H0)
# Summary statistics
summary(test_results$p_value)
# Histogram of p-values
hist(test_results$p_value, breaks = 20, col = "lightblue", main = "Histogram of p-values",
     xlab = "p-value", ylab = "Frequency", border = "black")

