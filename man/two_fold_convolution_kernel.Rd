% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_functions.R
\name{two_fold_convolution_kernel}
\alias{two_fold_convolution_kernel}
\title{Two-Fold Convolution Kernel Function}
\usage{
two_fold_convolution_kernel(u, kernel_func)
}
\arguments{
\item{u}{A numeric vector of points at which the two-fold convolution kernel is evaluated.}

\item{kernel_func}{A function representing the kernel to be convolved.}
}
\value{
A numeric vector of two-fold convolution kernel values corresponding to each input \code{u}.
}
\description{
This function computes the two-fold convolution of a given kernel function with itself.
The convolution is evaluated over a range of inputs \eqn{u} and is set to zero outside
the interval \([-2, 2]\).
}
\details{
The two-fold convolution kernel is defined as:
\deqn{
K^{(2)}(u) = \int_{-1}^{1} K(v) \cdot K(u - v) \, dv
}
where \eqn{K} is the original kernel function. The function evaluates this convolution for each
input \code{u} within the interval \([-2, 2]\) and sets it to zero outside this range.
}
\examples{
\dontrun{
# Define the Epanechnikov kernel function
epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}

# Define a range of u values
u_values <- seq(-3, 3, by = 0.1)

# Compute two-fold convolution kernel values
conv_kernel_values <- two_fold_convolution_kernel(u_values, kernel_func = epanechnikov_kernel)
}

}
\keyword{internal}
