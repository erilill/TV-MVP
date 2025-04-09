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
#' @section Authors and Maintainer:
#' Authors: Erik Lillrank and Yukai Yang \cr
#' Maintainer: Erik Lillrank \cr
#' Department of Statistics, Uppsala University \cr
#' \email{erik.lillrank.9607@@student.uu.se}, \email{yukai.yang@@statistik.uu.se}
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{determine_factors}}}{Selects the optimal number of factors via an information criterion.}
#'   \item{\code{\link{hyptest1}}}{Hypothesis test for time-invariant covariance matrices. Bootstrap p-values supported.}
#'   \item{\code{\link{rolling_time_varying_mvp}}}{Evaluates MVP performance in a rolling window framework.}
#'   \item{\code{\link{predict_portfolio}}}{Optimizes portfolio weights out-of-sample.}
#'   \item{\code{\link{time_varying_cov}}}{Estimates the time-varying covariance matrix.}
#' }
#'
#' @keywords package
"_PACKAGE"
