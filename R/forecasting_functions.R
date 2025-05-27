#' #' Rolling Window Time-Varying Minimum Variance Portfolio Optimization
#'
#' This function performs time-varying minimum variance portfolio (TV-MVP) optimization using
#' time-varying covariance estimation based on Local Principal Component Analysis (Local PCA). The
#' optimization is performed over a rolling window, with periodic rebalancing.
#' The procedure is available either as a stand-alone function or as a method in
#' the `TVMVP` R6 class.
#'
#' @param obj An object of class TVMVP with the data.
#' @param initial_window An integer specifying the number of periods used in the initial estimation window.
#' @param rebal_period An integer specifying the number of periods between portfolio rebalancing (e.g., weekly, monthly).
#' @param max_factors An integer indicating the maximum number of latent factors to consider in the factor model.
#' @param return_type A string indicating the return frequency. Options: `"daily"`, `"weekly"`, or `"monthly"`. Used for annualization of evaluation metrics.
#' @param kernel_func A kernel function to be used in the local PCA procedure. Default is `epanechnikov_kernel`.
#' @param rf A numeric vector or single value representing the risk-free rate. If `NULL`, it defaults to zero.
#' @param M0 An integer specifying the number of observations to leave out between the two sub-samples for the adaptive thresholding of the residual covariance estimation.
#' @param rho_grid A numeric sequence specifying the grid of rho values for threshold selection in covariance shrinkage. Default is `seq(0.005, 2, length.out = 30)`.
#' @param floor_value A small positive value to ensure numerical stability in the covariance matrix. Default is `1e-12`.
#' @param epsilon2 A small positive value used in the adaptive thresholding of the residual covariance estimation. Default is `1e-6`.
#'
#' @return An R6 object of class \code{RollingWindow} with the following accessible elements:
#' \describe{
#'   \item{\code{summary}}{A data frame of summary statistics for the TV-MVP and equal-weight portfolios, including cumulative excess return (CER), mean excess returns (MER), Sharpe ratio (SR), and standard deviation (SD) (raw and annualized).}
#'   \item{\code{TVMVP}}{A list containing rebalancing dates, estimated portfolio weights, and excess returns for the TV-MVP strategy.}
#'   \item{\code{Equal}}{A list with similar structure for the equal-weight portfolio.}
#' }
#'
#' @details
#' Two usage styles:
#' #' \preformatted{
#' # Function interface
#' results <- rolling_time_varying_mvp(
#'   obj = tv,
#'   initial_window = 50,
#'   rebal_period = 20,
#'   max_factors = 3,
#'   return_type = "daily",
#'   rf = NULL
#' )
#'
#' # R6 method interface
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#' results <- tv$rolling_time_varying_mvp(
#'   initial_window = 50,
#'   rebal_period = 20,
#'   max_factors = 3,
#'   return_type = "daily",
#'   rf = NULL)
#' }
#'
#' The function implements a rolling time-varying PCA approach to estimate latent factor structures
#' and uses a sparse residual covariance estimation method to improve covariance matrix estimation.
#' The covariance matrix is used to determine the global minimum variance portfolio (MVP), which is
#' rebalanced periodically according to the specified `rebal_period`. The number of factors is
#' determined by a BIC-type information criterion using the function `determine_factors`, updated
#' yearly. The bandwidth is determined by Silverman's rule of thumb, updated each rebalancing period.
#'
#' If `rf` is `NULL`, the risk-free rate is assumed to be zero.
#'
#' @section References: 
#' Lillrank, E. (2025). \ifelse{html}{
#'     \out{<a href='../doc/thesis.pdf'>A Time-Varying Factor Approach to Covariance Estimation</a>}
#'   }{Master’s thesis (PDF in inst/doc)}
#'   
#' Fan, Q., Wu, R., Yang, Y., & Zhong, W. (2024). Time-varying minimum variance portfolio. Journal of Econometrics, 239(2), 105339.
#' 
#' @examples
#' \donttest{
#' # Generate random returns for 20 assets over 100 periods
#' set.seed(123)
#' returns <- matrix(rnorm(20*100), nrow = 100, ncol = 20)
#' 
#' # Initialize object
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#'
#' # Run rolling TV-MVP optimization
#' results <- rolling_time_varying_mvp(
#'   obj = tv,
#'   initial_window = 50,
#'   rebal_period = 20,
#'   max_factors = 3,
#'   return_type = "daily",
#'   kernel_func = epanechnikov_kernel,
#'   rf = NULL
#' )
#'
#' # R6 method interface
#' results <- tv$rolling_time_varying_mvp(
#'   initial_window = 50,
#'   rebal_period = 20,
#'   max_factors = 3,
#'   return_type = "daily",
#'   rf = NULL)
#'
#' # Print summary
#' print(results)
#'
#' # Plot cumulative log excess returns
#' plot(results)
#' }
#' @importFrom stats var
#' @export
rolling_time_varying_mvp <- function(
    obj,
    initial_window,     # Number of periods in the initial estimation window
    rebal_period,       # Holding window length (HT in the paper)
    max_factors,
    return_type    = "daily",
    kernel_func    = epanechnikov_kernel,
    rf             = NULL,
    M0             = 10,                      # For covariance function
    rho_grid       = seq(0.005, 2, length.out = 30),  # For covariance function
    floor_value    = 1e-12,                   # For covariance function
    epsilon2       = 1e-6                     # For covariance function
) {
  if(inherits(obj, "TVMVP")){
    return(obj$rolling_time_varying_mvp(
      initial_window = initial_window, rebal_period = rebal_period,
      max_factors = max_factors, return_type    = return_type,
      kernel_func    = kernel_func, rf = rf,
      M0             = M0,
      rho_grid       = rho_grid,
      floor_value    = floor_value,
      epsilon2       = epsilon2
    ))
  } else {
    cli::cli_alert_danger("{.code obj} is not an object of {.strong TVMVP}")
  }
}

#' Predict Optimal Portfolio Weights Using Time-Varying Covariance Estimation
#'
#' This function estimates optimal portfolio weights using a time-varying covariance matrix
#' derived from Local Principal Component Analysis (Local PCA). The procedure is available either as a stand-alone
#' function or as a method in the `TVMVP` R6 class. It computes the following portfolios:
#' \enumerate{
#'   \item Global Minimum Variance Portfolio (MVP)
#'   \item Maximum Sharpe Ratio Portfolio (if \code{max_SR = TRUE})
#'   \item Return-Constrained Minimum Variance Portfolio (if \code{min_return} is provided)
#' }
#'
#' @param obj An object of class TVMVP with the data.
#' @param horizon Integer. Investment horizon over which expected return and risk are computed. Default is 1.
#' @param max_factors Integer. The number of latent factors to consider in the Local PCA model. Default is 3.
#' @param kernel_func Function. Kernel used for weighting observations in Local PCA. Default is \code{\link{epanechnikov_kernel}}.
#' @param min_return Optional numeric. If provided, the function computes a Return-Constrained Portfolio that targets this minimum return.
#' @param max_SR Logical. If TRUE, the Maximum Sharpe Ratio Portfolio is also computed. Default is \code{NULL}.
#' @param rf Numeric. Log risk-free rate. If \code{NULL}, defaults to 0.
#'
#' @return An object of class \code{PortfolioPredictions} (an R6 object) with:
#' \describe{
#'   \item{\code{summary}}{A data frame of evaluation metrics (expected return, risk, Sharpe ratio) for all computed portfolios.}
#'   \item{\code{MVP}}{A list containing the weights, expected return, risk, and Sharpe ratio for the Global Minimum Variance Portfolio.}
#'   \item{\code{max_SR}}{(Optional) A list with metrics for the Maximum Sharpe Ratio Portfolio.}
#'   \item{\code{MVPConstrained}}{(Optional) A list with metrics for the Return-Constrained Portfolio.}
#' }
#'
#' @section Methods:
#' The returned object includes:
#' \itemize{
#'   \item \code{$print()}: Nicely prints summary and portfolio access information.
#'   \item \code{$getWeights(method = c("MVP", "max_SR", "MVPConstrained"))}: Retrieves the weights for the selected portfolio.
#' }
#'
#' @details
#' Two usage styles:
#'
#' #' \preformatted{
#'
#' # R6 method interface
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#' tv$determine_factors(max_m=5)
#' prediction <- tv$predict_portfolio(horizon = 1, min_return = 0.01, max_SR = TRUE)
#' 
#' #' # Function interface
#' prediction <- predict_portfolio(obj, horizon = 5, m = 2, min_return = 0.01, max_SR=TRUE)
#' }
#' The methods can then be used on \code{prediction} to retrieve the weights.
#'
#' The function estimates a time-varying covariance matrix using Local PCA:
#' \deqn{\hat{\Sigma}_{r,t}=\hat{\Lambda}_t \hat{\Sigma}_F \hat{\Lambda}_t' + \tilde{\Sigma}_e}
#' Where \eqn{\hat{\Lambda}_t} is the factor loadings at time t, \eqn{\hat{\Sigma}_F} is the factor covariance matrix, and \eqn{\tilde{\Sigma}_e} is regularized covariance matrix of the idiosyncratic errors.
#'
#' It forecasts asset-level expected returns using a simple ARIMA model selection procedure and uses these in portfolio optimization.
#' Optimization strategies include:
#' \itemize{
#'   \item Global minimum variance (analytical)
#'   \item Maximum Sharpe ratio (if \code{max_SR = TRUE})
#'   \item Minimum variance with expected return constraint (if \code{min_return} is provided)
#' }
#'
#' @section References: 
#' Lillrank, E. (2025). \ifelse{html}{
#'     \out{<a href='../doc/thesis.pdf'>A Time-Varying Factor Approach to Covariance Estimation</a>}
#'   }{Master’s thesis (PDF in inst/doc)}
#'   
#' Fan, Q., Wu, R., Yang, Y., & Zhong, W. (2024). Time-varying minimum variance portfolio. Journal of Econometrics, 239(2), 105339.
#' 
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' returns <- matrix(rnorm(200 * 20, mean = 0, sd = 0.02), ncol = 20)
#' 
#' # Initialize object
#' tv <- TVMVP$new()
#' tv$set_data(returns)
#' 
#' # Optimize weights and predict returns
#' result <- predict_portfolio(
#'   tv,
#'   horizon = 5,
#'   m = 3,
#'   min_return = 0.02,
#'   max_SR = TRUE
#' )
#'
#' # Print the portfolio performance summary
#' print(result)
#'
#' # Access MVP weights
#' result$getWeights("MVP")
#'
#' # Access Max Sharpe weights (if computed)
#' result$getWeights("max_SR")
#'
#' # Or use R6 method interface
#' tv$determine_factors(max_m=5)
#' prediction <- tv$predict_portfolio(horizon = 1, min_return)
#' prediction
#' prediction$getWeights("MVPConstrained")
#' }
#' 
#' @export
predict_portfolio <- function(
    obj,
    horizon = 1,
    max_factors = 3,
    kernel_func = epanechnikov_kernel,
    min_return = NULL,
    max_SR = NULL,  # flag: if TRUE, compute maximum Sharpe portfolio
    rf = NULL
) {
  if(inherits(obj, "TVMVP")){
    return(obj$predict_portfolio(
      horizon = horizon, max_factors = max_factors, kernel_func = kernel_func,
      min_return = min_return, max_SR = max_SR, rf = rf))
  } else {
    cli::cli_alert_danger("{.code obj} is not an object of {.strong TVMVP}")
  }
}

#' Function to compute expected returns using a simple model selection approach
#' @param returns T times p matrix of returns
#' @param horizon Length of forecasting horizon
#' @importFrom stats arima AIC predict
comp_expected_returns <- function(returns, horizon) {
  exp_ret <- numeric(ncol(returns))

  for (i in seq_len(ncol(returns))) {
    candidate_models <- list()
    aics <- numeric()

    for (order in list(c(0,0,0), c(1,0,0), c(0,0,1), c(1,0,1))) {
      model <- tryCatch(
        arima(returns[, i], order = order),
        error = function(e) NULL
      )
      candidate_models <- c(candidate_models, list(model))
      aics <- c(aics, if (!is.null(model)) AIC(model) else Inf)
    }

    # If all models failed, fallback to mean return
    if (all(is.infinite(aics))) {
      exp_ret[i] <- mean(returns[, i])
    } else {
      best_model <- candidate_models[[which.min(aics)]]
      fc <- predict(best_model, n.ahead = horizon)$pred
      exp_ret[i] <- mean(fc)
    }
  }

  return(exp_ret)
}


#' @import R6
#' @import cli
PortfolioPredictions <- R6Class("PortfolioPredictions",
                                public = list(
                                  summary = NULL,
                                  MVP = NULL,
                                  max_SR = NULL,
                                  MVPConstrained = NULL,

                                  initialize = function(summary, MVP, max_SR = NULL, MVPConstrained = NULL) {
                                    self$summary <- summary
                                    self$MVP <- MVP
                                    self$max_SR <- max_SR
                                    self$MVPConstrained <- MVPConstrained
                                  },

                                  print = function(...) {
                                    cli::cli_h1("Portfolio Optimization Predictions")
                                    cli::cli_rule()
                                    cli::cli_h2("Summary Metrics")
                                    df <- self$summary
                                    df$Method <- with(df, ifelse(Method == "MVP",
                                                                 "Minimum Variance Portfolio",
                                                                 ifelse(Method == "max_SR",
                                                                        "Maximum SR Portfolio",
                                                                        ifelse(Method == "MVPConstrained",
                                                                               "Return-Constrained Portfolio", Method))))
                                    print(df, row.names = FALSE)
                                    cli::cli_rule()
                                    cli::cli_h2("Detailed Components")
                                    cli::cli_text("The detailed portfolio outputs are stored in the following elements:")
                                    cli::cli_text("  - MVP: Use object$MVP")
                                    if (!is.null(self$max_SR)) {
                                      cli::cli_text("  - Maximum Sharpe Ratio Portfolio: Use object$max_SR")
                                    }
                                    if (!is.null(self$MVPConstrained)) {
                                      cli::cli_text("  - Minimum Variance Portfolio with Return Constraint: Use object$MVPConstrained")
                                    }
                                    invisible(self)
                                  },

                                  getWeights = function(method = "MVP") {
                                    switch(method,
                                           MVP = self$MVP$weights,
                                           max_SR = {
                                             if (!is.null(self$max_SR)) {
                                               self$max_SR$weights
                                             } else {
                                               cli::cli_alert_danger("max_SR portfolio not available!")
                                               NULL
                                             }
                                           },
                                           MVPConstrained = {
                                             if (!is.null(self$MVPConstrained)) {
                                               self$MVPConstrained$weights
                                             } else {
                                               cli::cli_alert_danger("MVPConstrained portfolio not available!")
                                               NULL
                                             }
                                           },
                                           stop("Method not found")
                                    )
                                  }
                                )
)
#' @import R6
#' @import cli
RollingWindow <- R6::R6Class(
  "RollingWindow",
  public = list(
    summary = NULL,
    TVMVP   = NULL,
    Equal  = NULL,

    initialize = function(summary, TVMVP, Equal) {
      self$summary <- summary
      self$TVMVP   <- TVMVP
      self$Equal   <- Equal
    },

    print = function(...) {
      # Header
      cli::cli_h1("Rolling Window Portfolio Analysis")
      cli::cli_rule()

      # Print summary
      cli::cli_h2("Summary Metrics")
      df <- self$summary

      # Here, if you want, you can rename the 'Method' column so it prints
      # nicely, or leave it as is. For example:
      # df$Method <- with(df, ifelse(Method == "Time-Varying MVP",
      #                              "Time-Varying MVP Portfolio",
      #                       ifelse(Method == "Equal Weight",
      #                              "Equal-Weighted Portfolio",
      #                              Method)))

      print(df, row.names = FALSE)

      cli::cli_rule()
      cli::cli_h2("Detailed Components")
      cli::cli_text("The detailed portfolio outputs are stored in the following elements:")
      cli::cli_text("  - Time-Varying MVP: Access via `$TVMVP`")
      cli::cli_text("  - Equal Weight: Access via `$Equal`")

      invisible(self)
    },

    getWeights = function(method = c("TVMVP", "Equal")) {
      method <- match.arg(method)
      if (method == "TVMVP") {
        return(self$TVMVP$weights)
      } else {
        return(self$Equal$weights)
      }
    }
  )
)
#' @importFrom graphics legend lines par plot
#' @export
#' @method plot RollingWindow
plot.RollingWindow <- function(x, ...) {
  stopifnot(inherits(x, "RollingWindow"))

  # Calculate cumulative returns
  tvmvp_cum <- cumsum(x$TVMVP$returns)
  equal_cum <- cumsum(x$Equal$returns)
  T_len <- length(tvmvp_cum)

  # Adjust margins to make room for legend
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))  # Reset par on exit
  par(mar = c(6, 4, 4, 2) + 0.1)  # Extra space at bottom

  # Set up plotting area
  plot(
    1:T_len, tvmvp_cum,
    type = "l", col = "blue", lwd = 2,
    ylim = range(c(tvmvp_cum, equal_cum)),
    xlab = "Time",
    ylab = "Cumulative log Excess Return",
    main = "Cumulative log Excess Returns: TVMVP vs Equal",
    ...
  )

  # Add Equal Weight portfolio line
  lines(1:T_len, equal_cum, col = "red", lwd = 2, lty = 2)

  # Add legend BELOW the plot
  legend(
    "bottom", inset = -0.45, xpd = TRUE,
    legend = c("Time-Varying MVP", "Equal Weight"),
    col = c("blue", "red"),
    lty = c(1, 2), lwd = 2,
    horiz = TRUE, bty = "n", cex = 0.9
  )
}

