# R/evaluation_functions.R

#' Evaluate Portfolio Forecasts Using Various Error Metrics
#'
#' This function assesses the accuracy and reliability of portfolio forecasts by computing
#' multiple error metrics. It evaluates risk prediction accuracy, portfolio weight accuracy,
#' Sharpe ratio estimation accuracy, covariance matrix prediction accuracy, and residual
#' covariance matrix prediction accuracy.
#'
#' @param returns A numeric matrix of actual asset returns with \eqn{T} rows (time periods)
#'   and \eqn{p} columns (assets).
#' @param forecasts A numeric matrix of predicted means with \eqn{T} rows (time periods)
#'   and \eqn{p} columns (assets).
#' @param est_covariances A list of length \eqn{T}, where each element is a \eqn{p \times p}
#'   matrix of predicted covariances for each time period.
#' @param residual_covariances_pred A list of length \eqn{T}, where each element is a
#'   \eqn{p \times p} matrix of predicted residual covariances for each time period.
#' @param weights_est A list of length \eqn{T}, where each element is a numeric vector
#'   of portfolio weights of length \eqn{p} for each time period.
#' @param window_eval A numeric vector specifying the time indices to evaluate, e.g., \code{(W + 1):T}.
#' @param realized_covariances A list of length \eqn{T}, where each element is a
#'   \eqn{p \times p} matrix of realized covariances for each time period. Defaults to \code{NULL}.
#' @param realized_resid_covariances A list of length \eqn{T}, where each element is a
#'   \eqn{p \times p} matrix of realized residual covariances for each time period. Defaults to \code{NULL}.
#' @param true_weights A list of length \eqn{T}, where each element is a numeric vector
#'   of true portfolio weights of length \eqn{p} for each time period. Defaults to \code{NULL}.
#' @param realized_sharpes A numeric vector of length \eqn{T} representing the realized
#'   Sharpe Ratios for each time period. Defaults to \code{NULL}.
#'
#' @return A list containing the following error metrics:
#' \describe{
#'   \item{\code{risk_error}}{Numeric scalar. The average absolute difference between predicted
#'     and realized portfolio risks over the evaluation window.}
#'   \item{\code{weight_error}}{Numeric scalar. The average Euclidean distance between estimated
#'     and true portfolio weights over the evaluation window.}
#'   \item{\code{sharpe_error}}{Numeric scalar. The average absolute difference between
#'     predicted and realized Sharpe Ratios over the evaluation window.}
#'   \item{\code{covariance_error}}{Numeric scalar. The average Frobenius norm of the difference
#'     between predicted and realized covariance matrices over the evaluation window.}
#'   \item{\code{residual_cov_error}}{Numeric scalar. The average Frobenius norm of the difference
#'     between predicted and realized residual covariance matrices over the evaluation window.}
#' }
#'
#' @details
#' The function computes the following error metrics over the specified evaluation window:
#' \enumerate{
#'   \item \strong{Risk Error:} For each portfolio, it calculates the absolute difference between the
#'     predicted risk (standard deviation) and the realized risk, then averages these differences.
#'   \item \strong{Weight Error:} For each portfolio, it computes the Euclidean distance between the
#'     estimated portfolio weights and the true portfolio weights, then averages these distances.
#'   \item \strong{Sharpe Ratio Error:} For each portfolio, it calculates the absolute difference
#'     between the predicted Sharpe Ratio (based on forecasted means and estimated covariances)
#'     and the realized Sharpe Ratio, then averages these differences.
#'   \item \strong{Covariance Error:} For each covariance matrix, it computes the Frobenius norm of
#'     the difference between the predicted and realized covariance matrices, then averages these norms.
#'   \item \strong{Residual Covariance Error:} Similarly, for each residual covariance matrix, it
#'     computes the Frobenius norm of the difference between the predicted and realized residual
#'     covariance matrices, then averages these norms.
#' }
#'
#' @examples
#' # Example parameters
#' T <- 100  # Number of time periods
#' p <- 50   # Number of assets
#'
#' # Simulate actual returns
#' set.seed(123)
#' returns <- matrix(rnorm(T * p, mean=0.001, sd=0.02), nrow=T, ncol=p)
#'
#' # Simulate forecasts (predicted means)
#' forecasts <- matrix(rnorm(T * p, mean=0.001, sd=0.015), nrow=T, ncol=p)
#'
#' # Simulate estimated covariances
#' est_covariances <- lapply(1:T, function(t) {
#'   diag(runif(p, min=0.0001, max=0.002)) + matrix(runif(p^2, min=-0.0005, max=0.0005), nrow=p)
#' })
#'
#' # Simulate predicted residual covariances
#' residual_covariances_pred <- lapply(1:T, function(t) {
#'   diag(runif(p, min=0.0001, max=0.002))
#' })
#'
#' # Simulate estimated weights
#' weights_est <- lapply(1:T, function(t) {
#'   w <- runif(p)
#'   w / sum(w)
#' })
#'
#' # Define evaluation window
#' window_eval <- 51:100
#'
#' # Simulate realized covariances
#' realized_covariances <- lapply(1:T, function(t) {
#'   cov(returns[(max(1, t-10)):t, , drop=FALSE])
#' })
#'
#' # Simulate realized residual covariances
#' realized_resid_covariances <- lapply(1:T, function(t) {
#'   diag(runif(p, min=0.0001, max=0.002))
#' })
#'
#' # Simulate true weights
#' true_weights <- lapply(1:T, function(t) {
#'   w <- runif(p)
#'   w / sum(w)
#' })
#'
#' # Simulate realized Sharpe Ratios
#' realized_sharpes <- runif(T, min=0.5, max=2.0)
#'
#' # Evaluate forecasts
#' metrics <- evaluate_forecasts(
#'   returns = returns,
#'   forecasts = forecasts,
#'   est_covariances = est_covariances,
#'   residual_covariances_pred = residual_covariances_pred,
#'   weights_est = weights_est,
#'   window_eval = window_eval,
#'   realized_covariances = realized_covariances,
#'   realized_resid_covariances = realized_resid_covariances,
#'   true_weights = true_weights,
#'   realized_sharpes = realized_sharpes
#' )
#' print(metrics)
#'
#' @export
evaluate_forecasts <- function(
    returns,
    forecasts,                  # T x p matrix of predicted means
    est_covariances,            # list of length T, each a p x p predicted covariance
    residual_covariances_pred,  # list of length T, each a p x p predicted residual covariance
    weights_est,                # list of length T, each a p-vector of portfolio weights
    window_eval,                # Vector of time indices, e.g., (W +1):T
    realized_covariances   = NULL,       # list of length T
    realized_resid_covariances = NULL,   # list of length T
    true_weights           = NULL,       # list of length T
    realized_sharpes       = NULL        # numeric vector of length T
) {
  # Initialize metrics
  metrics <- list(
    risk_error = NA,
    weight_error = NA,
    sharpe_error = NA,
    covariance_error = NA,
    residual_cov_error = NA
  )

  # 1. Risk Error
  if (!is.null(realized_covariances) && !is.null(weights_est)) {
    risk_diffs <- mapply(function(w, est_cov, real_cov) {
      if (is.null(w) || is.null(est_cov) || is.null(real_cov)) return(NA)
      pred_risk <- sqrt(as.numeric(t(w) %*% est_cov %*% w))
      real_risk <- sqrt(as.numeric(t(w) %*% real_cov %*% w))
      return(abs(real_risk - pred_risk))
    },
    w = weights_est[window_eval],
    est_cov = est_covariances[window_eval],
    real_cov = realized_covariances[window_eval])

    metrics$risk_error <- mean(risk_diffs, na.rm = TRUE)
  }

  # 2. Weight Error
  if (!is.null(true_weights) && !is.null(weights_est)) {
    weight_diffs <- mapply(function(w_est, w_true) {
      if (is.null(w_est) || is.null(w_true)) return(NA)
      return(sqrt(sum((w_est - w_true)^2)))
    },
    w_est = weights_est[window_eval],
    w_true = true_weights[window_eval])

    metrics$weight_error <- mean(weight_diffs, na.rm = TRUE)
  }

  # 3. Sharpe Ratio Error
  if (!is.null(realized_sharpes) && !is.null(weights_est) && !is.null(forecasts)) {
    sharpe_diffs <- mapply(function(w, forecast_row, realized_sr, est_cov) {
      if (is.null(w) || is.null(forecast_row) || is.na(realized_sr) || is.null(est_cov)) return(NA)
      mu_hat <- sum(forecast_row * w, na.rm = TRUE)
      sigma_hat <- sqrt(as.numeric(t(w) %*% est_cov %*% w))
      if (sigma_hat == 0 || is.na(sigma_hat)) return(NA)
      sr_hat <- mu_hat / sigma_hat
      return(abs(realized_sr - sr_hat))
    },
    w = weights_est[window_eval],
    forecast_row = split(forecasts[window_eval, ], row(forecasts[window_eval, ])),
    realized_sr = realized_sharpes[window_eval],
    est_cov = est_covariances[window_eval])

    metrics$sharpe_error <- mean(sharpe_diffs, na.rm = TRUE)
  }

  # 4. Covariance Error
  if (!is.null(realized_covariances) && !is.null(est_covariances)) {
    cov_diffs <- mapply(function(est_cov, real_cov) {
      if (is.null(est_cov) || is.null(real_cov)) return(NA)
      return(norm(est_cov - real_cov, type = "F"))
    },
    est_cov = est_covariances[window_eval],
    real_cov = realized_covariances[window_eval])

    metrics$covariance_error <- mean(cov_diffs, na.rm = TRUE)
  }

  # 5. Residual Covariance Error
  if (!is.null(realized_resid_covariances) && !is.null(residual_covariances_pred)) {
    resid_cov_diffs <- mapply(function(pred_resid, real_resid) {
      if (is.null(pred_resid) || is.null(real_resid)) return(NA)
      return(norm(pred_resid - real_resid, type = "F"))
    },
    pred_resid = residual_covariances_pred[window_eval],
    real_resid = realized_resid_covariances[window_eval])

    metrics$residual_cov_error <- mean(resid_cov_diffs, na.rm = TRUE)
  }

  return(metrics)
}
