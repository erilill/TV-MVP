#' Epanechnikov Kernel Function
#'
#' This function computes the value of the Epanechnikov kernel for a given input \eqn{u}.
#' The Epanechnikov kernel is a popular choice in kernel density estimation due to its optimal
#' properties in minimizing mean integrated squared error.
#'
#' @param u A numeric vector of points at which the kernel is evaluated.
#'
#' @return A numeric vector of kernel values corresponding to each input \code{u}.
#'
#' @details
#' The Epanechnikov kernel is defined as:
#' \deqn{
#' K(u) = \begin{cases}
#' \frac{3}{4}(1 - u^2) & \text{if } |u| \leq 1, \\
#' 0 & \text{otherwise}.
#' \end{cases}
#' }
#'
#' This function applies the above formula to each element of the input vector \code{u}.
#'
#' @examples
#' \dontrun{
#' # Define a range of u values
#' u_values <- seq(-1.5, 1.5, by = 0.1)
#'
#' # Compute Epanechnikov kernel values
#' kernel_values <- epanechnikov_kernel(u_values)
#'
#' }
#' @export
epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}
#' Boundary Kernel Function
#'
#' This function computes boundary kernel weights for a given time period \eqn{t} within a dataset
#' of size \eqn{T}. It adjusts the kernel weights near the boundaries to account for edge effects,
#' ensuring that the weights sum to one.
#'
#' @param t An integer specifying the current time period for which the kernel weights are computed.
#' @param r An integer representing the reference time period.
#' @param T An integer indicating the total number of time periods in the dataset.
#' @param h A numeric value representing the bandwidth parameter for the kernel function.
#' @param kernel_func A function representing the kernel used for weighting.
#'
#' @return A numeric scalar representing the boundary-adjusted kernel weight for the given time period.
#'
#' @details
#' The boundary kernel function adjusts kernel weights near the start and end of the dataset to mitigate
#' edge effects commonly encountered in kernel-based methods. The function performs the following steps:
#' \enumerate{
#'   \item Scales the difference between the current time \code{t} and reference time \code{r} by the
#'   product of total time periods \code{T} and bandwidth \code{h}.
#'   \item Applies the kernel function to the scaled difference and adjusts by the bandwidth.
#'   \item Determines if the current time period is within the lower or upper boundary regions based on
#'   \eqn{T_h = \lfloor T \times h \rfloor}.
#'   \item Computes the integral of the kernel function over the adjusted limits to ensure the weights
#'   sum to one in boundary regions.
#' }
#'
#' @examples
#' \dontrun{
#' # Define the Epanechnikov kernel function
#' epanechnikov_kernel <- function(u) {
#'   ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
#' }
#'
#' # Compute boundary kernel weight for time period t = 5, reference r = 5, T = 100, h = 0.1
#' weight <- boundary_kernel(t = 5, r = 5, T = 100, h = 0.1, kernel_func = epanechnikov_kernel)
#' print(weight)
#' }
#' @export
boundary_kernel <- function(t, r, iT, h, kernel_func) {
  scaled_diff <- (t - r) / (iT * h)
  k_val <- kernel_func(scaled_diff) / h

  # Determine the region of r
  Th_floor <- floor(iT * h)

  if (r < Th_floor) {
    # Lower boundary case
    integral_val <- integrate(kernel_func, lower = -r / (iT * h), upper = 1)$value
    return(k_val / integral_val)
  } else if (r > (iT - Th_floor)) {
    # Upper boundary case
    integral_val <- integrate(kernel_func, lower = -1, upper = (1 - r / iT) / h)$value
    return(k_val / integral_val)
  } else {
    # Middle region
    return(k_val)
  }
}
#' Two-Fold Convolution Kernel Function
#'
#' This function computes the two-fold convolution of a given kernel function with itself.
#' The convolution is evaluated over a range of inputs \eqn{u} and is set to zero outside
#' the interval \([-2, 2]\).
#'
#' @param u A numeric vector of points at which the two-fold convolution kernel is evaluated.
#' @param kernel_func A function representing the kernel to be convolved. 
#'
#' @return A numeric vector of two-fold convolution kernel values corresponding to each input \code{u}.
#'
#' @details
#' The two-fold convolution kernel is defined as:
#' \deqn{
#' K^{(2)}(u) = \int_{-1}^{1} K(v) \cdot K(u - v) \, dv
#' }
#' where \eqn{K} is the original kernel function. The function evaluates this convolution for each
#' input \code{u} within the interval \([-2, 2]\) and sets it to zero outside this range.
#'
#' @examples
#' \dontrun
#' # Define the Epanechnikov kernel function
#' epanechnikov_kernel <- function(u) {
#'   ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
#' }
#'
#' # Define a range of u values
#' u_values <- seq(-3, 3, by = 0.1)
#'
#' # Compute two-fold convolution kernel values
#' conv_kernel_values <- two_fold_convolution_kernel(u_values, kernel_func = epanechnikov_kernel)
#'
#' }
#'
#' @export
two_fold_convolution_kernel <- function(u, kernel_func) {
  result <- ifelse(
    abs(u) <= 2,
    sapply(u, function(u_val) {
      integrand <- function(v) kernel_func(v) * kernel_func(u_val - v)
      integrate(integrand, lower = -1, upper = 1)$value
    }),
    0
  )
  return(result)
}

