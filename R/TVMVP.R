################################################################################
## Package name: TVMVP
## Authors: Erik Lillrank and Yukai Yang
## Department of Statistics, Uppsala University
## April 2025
################################################################################

#' TVMVP: Time-Varying Minimum Variance Portfolio Optimization
#'
#' The TVMVP package provides tools for estimating time-dependent covariance 
#' matrices using kernel-weighted principal component analysis. These estimates 
#' can then be used for portfolio optimization in a dynamic setting.
#'
#' The method involves five steps: (1) determining the number of factors, 
#' (2) estimating kernel-weighted PCA, (3) regularizing the idiosyncratic error 
#' covariance, (4) estimating the total covariance matrix, and (5) computing 
#' optimal portfolio weights.
#'
#' An optional step includes a hypothesis test to check whether the covariance 
#' matrix is time-invariant.
#' 
#'
#' @section Authors and Maintainer:
#' Authors: Erik Lillrank and Yukai Yang \cr
#' Maintainer: Erik Lillrank \cr
#' Department of Statistics, Uppsala University \cr
#' \email{erik.lillrank@@gmail.com}, \email{yukai.yang@@statistik.uu.se}
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{determine_factors}}}{Selects the optimal number of factors via an information criterion.}
#'   \item{\code{\link{hyptest1}}}{Hypothesis test for time-invariant covariance matrices. Bootstrap p-values supported.}
#'   \item{\code{\link{predict_portfolio}}}{Optimizes portfolio weights for out-of-sample prediction of portfolio performance.}
#'   \item{\code{\link{rolling_time_varying_mvp}}}{Evaluates MVP performance in a rolling window framework.}
#'   \item{\code{\link{time_varying_cov}}}{Estimates the time-varying covariance matrix.}
#'   \item{\code{\link{silverman}}}{Silverman's rule of thumb bandwidth formula.}
#' }
#' 
#' 
#' @section Class:
#' \describe{
#'  \item{\code{\link{TVMVP}}}{Time Varying Minimum Variance Portfolio (TVMVP) Class.}
#' } 
#' 
#'
#' @details
#' The default kernel function in the package is the Epanechnikov kernel. Other
#' kernel functions can also be used, however these are not implemented in the
#' package. In order to do this, write an R function with an integrable
#' kernel function and use this as input in the functions with argument
#' \code{kernel_func}. It should be constructed as \code{custom_kernel <- function(u){...}}.
#' 
#' Similarly, the bandwidth function which is implemented in the package is 
#' the Silverman's rule of thumb. For most functions, simply set \code{bandwidth}
#' to your preferred bandwidth, however for \code{\link{rolling_time_varying_mvp}},
#' only Silverman's is implemented in this version of the package.
#' 
#' @section References: 
#' Lillank, E. (2025). A Time-Varying Factor Approach to Covariance Estimation.
#' 
#' @keywords package
"_PACKAGE"