test_that("predict_portfolio returns a PortfolioPredictions R6 object", {
  set.seed(123)
  T <- 100  # Number of time periods
  p <- 20   # Number of assets
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

  tmp <- TVMVP$new()
  tmp$set(data = returns)

  prediction <- predict_portfolio(
    tmp,
    horizon = 21,
    kernel_func = epanechnikov_kernel,
    max_factors = 10,
    min_return = 0.5
  )

  # R6 class check
  expect_true(inherits(prediction, "PortfolioPredictions"))

  # Optional: check for something else
})
