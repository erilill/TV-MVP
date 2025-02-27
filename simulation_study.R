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


library(quantmod) # API to get data from yahoo finance

# Tickers of a handull of stocks with no missing values between 2022-06-30 - 2023-06-30
# These where the tickers I could easily get, it is about 40% of all stocks on stockholmsbörsen.
# There was no cherry-picking so I think it is fine. 
swe_stocks <- c("ATCO-A.ST",
                "ATCO-B.ST",
                "SCA-A.ST",
                "SCA-B.ST",
                "SKA-B.ST",
                "ELUX-A.ST",
                "ELUX-B.ST",
                "HM-B.ST",
                "HUFV-A.ST",
                "ERIC-A.ST",
                "ERIC-B.ST",
                "HOLM-A.ST",
                "HOLM-B.ST",
                "SAND.ST",
                "SKF-A.ST",
                "SKF-B.ST",
                "TREL-B.ST",
                "VOLV-A.ST",
                "BERG-B.ST",
                "HEXA-B.ST",
                "LATO-B.ST",
                "INDU-A.ST",
                "INDU-C.ST",
                "INVE-A.ST",
                "INVE-B.ST",
                "RATO-A.ST",
                "RATO-B.ST",
                "ORES.ST",
                "JM.ST",
                "NCC-A.ST",
                "NCC-B.ST",
                "AZA.ST",
                "CAT-A.ST",
                "CAT-B.ST",
                "SINT.ST",
                "NOKIA-SEK.ST",
                "GETI-B.ST",
                "SEB-A.ST",
                "SEB-C.ST",
                "SHB-A.ST",
                "SHB-B.ST",
                "LUND-B.ST",
                "BEIJ-B.ST",
                "OEM-B.ST",
                "KABE-B.ST",
                "BILI-A.ST",
                "NOLA-B.ST",
                "WALL-B.ST",
                "SVOL-A.ST",
                "SVOL-B.ST",
                "BURE.ST",
                "AFRY.ST",
                "ACTI.ST",
                "PEAB-B.ST",
                "KINV-A.ST",
                "KINV-B.ST",
                "STWK.ST",
                "BEIA-B.ST",
                "FABG.ST",
                "VBG-B.ST",
                "RROS.ST",
                "BONG.ST",
                "HAV-B.ST",
                "ELAN-B.ST",
                "SSAB-A.ST",
                "SSAB-B.ST",
                "ENEA.ST",
                "RO.ST",
                "RAY-B.ST",
                "SECU-B.ST",
                "VSSAB-B.ST",
                "DORO.ST",
                "FPAR-A.ST",
                "EKTA-B.ST",
                "CNCJO-B.ST",
                "HEBA-B.ST",
                "SKIS-B.ST",
                "ATRLJ-B.ST",
                "ASSA-B.ST",
                "PRIC-B.ST",
                "ORTI-A.ST",
                "ORTI-B.ST",
                "MVIR.ST",
                "BETS-B.ST",
                "TEL2-A.ST",
                "TEL2-B.ST",
                "DURC-B.ST",
                "ALIV-SDB.ST",
                "FAG.ST",
                "CAST.ST",
                "NIBE-B.ST",
                "PROF-B.ST",
                "CORE-A.ST",
                "LAMM-B.ST",
                "AXFO.ST",
                "TRAC-B.ST",
                "MTG-A.ST",
                "MTG-B.ST",
                "SVED-B.ST",
                "PACT.ST",
                "KNOW.ST",
                "CTT.ST",
                "NDA-SE.ST",
                "NEWA-B.ST",
                "JLT.ST",
                "ADDV-A.ST",
                "FING-B.ST",
                "SAAB-B.ST",
                "BIOG-B.ST",
                "PREV-B.ST",
                "SWEC-A.ST",
                "SWEC-B.ST",
                "SOF-B.ST",
                "PROB.ST",
                "STE-A.ST",
                "STE-R.ST",
                "IAR-B.ST",
                "SECT-B.ST",
                "ARBO-A.ST",
                "AZN.ST",
                "SAGA-A.ST",
                "IS.ST",
                "NETI-B.ST",
                "ANOD-B.ST",
                "MSON-A.ST",
                "MSON-B.ST",
                "ABB.ST",
                "PION-B.ST",
                "VIT-B.ST",
                "NTEK-B.ST",
                "TIETOS.ST",
                "CLAS-B.ST",
                "BALD-B.ST",
                "MSAB-B.ST",
                "PREC.ST",
                "HIFA-B.ST",
                "NGS.ST",
                "STRAX.ST",
                "MYCR.ST",
                "ANOT.ST",
                "MEKO.ST",
                "EPEN.ST",
                "TELIA.ST",
                "MOMENT.ST",
                "BIOT.ST",
                "STAR-B.ST",
                "ENRO.ST",
                "SVIK.ST",
                "BTS-B.ST",
                "AQ.ST",
                "BINV.ST",
                "VITR.ST",
                "COALA.ST",
                "LAGR-B.ST",
                "ADDT-B.ST",
                "ORRON.ST",
                "BILL.ST",
                "BOL.ST",
                "ALFA.ST",
                "INTRUM.ST",
                "NOBI.ST",
                "REJL-B.ST",
                "LUMI.ST",
                "TIGO-SDB.ST",
                "TETY.ST",
                "ITAB.ST",
                "IVSO.ST",
                "NOTE.ST",
                "MCAP.ST",
                "BORG.ST",
                "FPIP.ST",
                "WIHL.ST",
                "AAK.ST",
                "INDT.ST",
                "TRAD.ST",
                "ORX.ST",
                "VIVE.ST",
                "CEVI.ST",
                "CATE.ST",
                "ELON.ST",
                "DIOS.ST",
                "HUSQ-B.ST",
                "HUSQ-A.ST",
                "SOBI.ST",
                "G5EN.ST",
                "BEGR.ST",
                "LIAB.ST",
                "WISE.ST",
                "NAXS.ST",
                "FNOX.ST",
                "NMAN.ST",
                "LOGI-A.ST",
                "VNV.ST",
                "CRAD-B.ST",
                "RAIL.ST",
                "SYSR.ST",
                "HNSA.ST",
                "HMS.ST",
                "EAST.ST",
                "DUNI.ST",
                "HPOL-B.ST",
                "HTRO.ST",
                "SEZI.ST",
                "VESTUM.ST",
                "EWRK.ST",
                "ABLI.ST",
                "CLA-B.ST",
                "LOOMIS.ST",
                "EOLU-B.ST",
                "CORE-PREF.ST",
                "ARISE.ST",
                "BMAX.ST",
                "AOI.ST",
                "EPIS-B.ST",
                "NELLY.ST",
                "KDEV.ST",
                "DEDI.ST",
                "BULTEN.ST",
                "MOB.ST",
                "BOUL.ST",
                "STEF-B.ST",
                "SIVE.ST",
                "CCC.ST",
                "LUC.ST",
                "CRED-A.ST",
                "FASTAT.ST",
                "MANG.ST",
                "XVIVO.ST",
                "EGTX.ST",
                "ARP.ST",
                "SAGA-B.ST",
                "IMMU.ST",
                "PLAZ-B.ST",
                "BUFAB.ST",
                "SANION.ST",
                "HANZA.ST",
                "SCST.ST",
                "BACTI-B.ST",
                "FOI-B.ST",
                "INWI.ST",
                "GRNG.ST",
                "SBB-B.ST",
                "LIFCO-B.ST",
                "THULE.ST",
                "NP3.ST",
                "PCELL.ST",
                "LUG.ST",
                "ELTEL.ST",
                "DUST.ST",
                "ACRI-A.ST",
                "SDIP-PREF.ST",
                "CANTA.ST",
                "EVO.ST",
                "HOFI.ST",
                "TROAX.ST",
                "K2A-PREF.ST",
                "TOBII.ST",
                "TRANS.ST",
                "VOLO-PREF.ST",
                "COOR.ST",
                "NIL-B.ST",
                "STAR-A.ST",
                "ALIG.ST",
                "PNDX-B.ST",
                "VEFAB.ST",
                "SINCH.ST",
                "BRAV.ST",
                "NICA.ST",
                "DOM.ST",
                "CAMX.ST",
                "ATT.ST",
                "SHOT.ST",
                "IMMNOV.ST",
                "SF.ST",
                "VICO.ST",
                "XBRANE.ST",
                "CTM.ST",
                "GARO.ST",
                "ALIF-B.ST",
                "HUM.ST",
                "IBT-B.ST",
                "RESURS.ST",
                "NWG.ST",
                "BONAV-A.ST",
                "BONAV-B.ST",
                "B3.ST",
                "TFBANK.ST",
                "ACAD.ST",
                "SYNACT.ST",
                "MAHA-A.ST",
                "BRIN-B.ST",
                "BICO.ST",
                "SAGA-D.ST",
                "EMBRAC-B.ST",
                "ATORX.ST",
                "VOLO.ST",
                "ONCO.ST",
                "IRLAB-A.ST",
                "MIPS.ST",
                "AMBEA.ST",
                "ISOFOL.ST",
                "ATIC.ST",
                "FMM-B.ST",
                "IPCO.ST",
                "INSTAL.ST",
                "SDIP-B.ST",
                "MTRS.ST",
                "MCOV-B.ST",
                "BOOZT.ST",
                "ESSITY-A.ST",
                "ESSITY-B.ST",
                "ALLIGO-B.ST",
                "BONEX.ST",
                "SEDANA.ST",
                "TRIAN-B.ST",
                "XSPRAY.ST",
                "BALCO.ST",
                "BIOA-B.ST",
                "FNM.ST",
                "SEAF.ST",
                "ARJO-B.ST",
                "CORE-B.ST",
                "NP3-PREF.ST",
                "CIBUS.ST",
                "GREEN.ST",
                "BHG.ST",
                "INFREA.ST",
                "OVZON.ST",
                "NCAB.ST",
                "BETCO.ST",
                "ARION-SDB.ST",
                "EPI-A.ST",
                "EPI-B.ST",
                "PENG-B.ST",
                "NYF.ST",
                "LIME.ST",
                "QLINEA.ST",
                "SBB-D.ST",
                "ACE.ST",
                "VPLAY-A.ST",
                "VPLAY-B.ST",
                "KAR.ST",
                "JOMA.ST",
                "K2A-B.ST",
                "8TRA.ST",
                "EQT.ST",
                "KFAST-B.ST",
                "FPAR-D.ST",
                "EPRO-B.ST",
                "GPG.ST",
                "READ.ST",
                "QLIRO.ST",
                "WBGR-B.ST",
                "NPAPER.ST",
                "SAVE.ST",
                "ANNE-B.ST",
                "FG.ST",
                "CINT.ST",
                "CS.ST",
                "PIERCE.ST",
                "ACRI-B.ST",
                "HEM.ST",
                "ARPL.ST",
                "MILDEF.ST",
                "LINC.ST",
                "SLEEP.ST",
                "RVRC.ST",
                "CORE-D.ST",
                "PRFO.ST",
                "CTEK.ST",
                "EMIL-PREF.ST",
                "STOR-B.ST",
                "TRUE-B.ST",
                "NETEL.ST",
                "NORB-B.ST",
                "VOLCAR-B.ST",
                "SYNSAM.ST",
                "LOGI-B.ST",
                "SFAB.ST",
                "KLARA-B.ST",
                "NIVI-B.ST",
                "NORVA.ST",
                "SLP-B.ST",
                "MMGR-B.ST",
                "EMIL-B.ST",
                "ENGCON-B.ST"
)


stock_env <- new.env() # avoid cluttering environment
getSymbols(swe_stocks, src = "yahoo", from = "2022-06-30", to = "2023-06-30", env = stock_env)
swe_stock_prices <- do.call(merge, lapply(swe_stocks, function(x) Cl(get(x, envir = stock_env))))
returns <- as.matrix(diff(log(swe_stock_prices))[-1])

pred <-predict_portfolio(returns[,1:100], 21, bandwidth_func = silverman, min_return = 0.05)
rolpred <- rolling_time_varying_mvp(returns[,1:100], 200, 5, 5)

## simulation using stock data

# Run PCA to get factors and loadings which will be used for simulation
swe_pca <- localPCA(returns[,1:100], silverman(returns[,1:100]), 2)

swe_sim <- function(swe_pca){
  for (t in nrow(swe_pca$factors[[1]])){
    residual_cov <- 
    pred <- tcrossprod(swe_pca$factors[[t]], swe_pca$loadings)
    
  }
}
