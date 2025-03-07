#' Compute Bandwidth Parameter Using Silverman's Rule of Thumb
#'
#' This function calculates the bandwidth parameter for kernel functions using Silverman's rule of thumb,
#' which is commonly used in kernel density estimation to determine an appropriate bandwidth.
#'
#' @param returns A numeric matrix of asset returns with \eqn{T} rows (time periods) and \eqn{p} columns (assets).
#'
#' @return A numeric value representing the computed bandwidth parameter based on Silverman's rule.
#'
#' @details
#' Silverman's rule of thumb for bandwidth selection is given by:
#' \deqn{bandwidth = \frac{2.35}{\sqrt{12}} \times T^{-0.2} \times p^{-0.1}}
#' where \eqn{T} is the number of time periods and \eqn{p} is the number of assets.
#'
#' If the number of time periods \code{T} and the number of assets \code{p} are not provided,
#' the function extracts these values from the \code{returns} matrix.
#'
#' @examples
#' # Simulate data for 50 assets over 200 time periods
#' set.seed(123)
#' T <- 200
#' p <- 50
#' returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)
#'
#' # Compute bandwidth using Silverman's rule of thumb
#' bw <- silverman(returns)
#' print(bw)
#'
#' @export
silverman <- function(returns){
    ip <- ncol(returns)
    iT <- nrow(returns)
  bandwidth <- (2.35/sqrt(12)) * iT^(-0.2) * ip^(-0.1)
  return(bandwidth)
}
