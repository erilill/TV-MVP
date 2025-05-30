% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localPOET.R
\name{estimate_residual_cov_poet_local}
\alias{estimate_residual_cov_poet_local}
\title{Estimate Local Covariance}
\usage{
estimate_residual_cov_poet_local(
  localPCA_results,
  returns,
  M0 = 10,
  rho_grid = seq(0.005, 2, length.out = 30),
  floor_value = 1e-12,
  epsilon2 = 1e-06
)
}
\arguments{
\item{localPCA_results}{A list containing the results from local PCA, with components:
\itemize{
  \item \code{loadings}: a list where each element is a \eqn{p × m} matrix of factor loadings.
  \item \code{f_hat}: a \eqn{T × m} matrix of estimated factors.
  \item \code{weights}: a list of kernel weight vectors.
}}

\item{returns}{A numeric matrix of asset returns with dimensions \eqn{T × p}.}

\item{M0}{Integer. The number of observations to leave out between the two sub-samples in the adaptive thresholding procedure. Default is 10.}

\item{rho_grid}{A numeric vector of candidate shrinkage parameters \eqn{\rho} used in \code{adaptive_poet_rho()}. Default is \code{seq(0.005, 2, length.out = 30)}.}

\item{floor_value}{A small positive number specifying the lower bound for eigenvalues in the final positive semidefinite repair. Default is \code{1e-12}.}

\item{epsilon2}{A small positive tuning parameter for the adaptive thresholding. Default is \code{1e-6}.}
}
\value{
A list containing:
  \itemize{
    \item \code{best_rho}: The selected shrinkage parameter \eqn{\hat{\rho}_t} for the local residual covariance.
    \item \code{residual_cov}: The shrunk residual covariance matrix \eqn{\hat{\Sigma}_e(T)}.
    \item \code{total_cov}: The final estimated time-varying covariance matrix \eqn{\Sigma_R(t)}.
    \item \code{loadings}: The local factor loadings \eqn{\Lambda_t} from the local PCA.
    \item \code{naive_resid_cov}: The raw (unshrunk) residual covariance matrix.
  }
}
\description{
This internal function computes a time-varying covariance matrix estimate for a given
window of asset returns by combining factor-based and sparse residual covariance estimation.
It uses results from a local PCA to form residuals and then applies an adaptive thresholding
procedure (via \code{adaptive_poet_rho()}) to shrink the residual covariance.
}
\details{
The function follows these steps:

\enumerate{
  \item **Local Residuals:**
        Extract the local loadings \eqn{\Lambda_t} from the last element of \code{localPCA_results\$loadings} and
        factors \eqn{\hat{F}} from \code{localPCA_results\$f_hat}. Let \eqn{w_t} denote the corresponding kernel weights.
        The local residuals are computed as:
        \deqn{U_{\text{local}} = R - F \Lambda_t,}
        where \eqn{R} is the returns matrix.

  \item **Adaptive Thresholding:**
        The function calls \code{adaptive_poet_rho()} on \eqn{U_{\text{local}} }to select an optimal shrinkage parameter
        \eqn{\hat{\rho}_t}.

  \item **Residual Covariance Estimation:**
        The raw residual covariance is computed as:
        \deqn{S_{u,\text{raw}} = \frac{1}{T} U_{\text{local}}^\top U_{\text{local}},}
        and a threshold is set as:
        \deqn{\text{threshold} = \hat{\rho}_t × \text{mean}(|S_{u,\text{raw}}|),}
        where the mean is taken over the off-diagonal elements.
        Soft-thresholding is then applied to obtain the shrunk residual covariance matrix \eqn{\hat{S}_u}.

  \item **Total Covariance Estimation:**
        The final covariance matrix is constructed by combining the factor component with the shrunk residual covariance:
        \deqn{\Sigma_R(t) = \Lambda_t \left(\frac{F^\top F}{T}\right) \Lambda_t^\top + \hat{S}_u.}

  \item **PSD Repair:**
        A final positive semidefinite repair is performed by flooring eigenvalues at \code{floor_value} and symmetrizing the matrix.
}
}
\keyword{internal}
