test_that("expanding_tvmvp returns a ExpandingWindow R6 object", {
  set.seed(123)
  T <- 100  # Number of time periods
  p <- 20   # Number of assets
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

  tmp <- TVMVP$new()
  tmp$set(data = returns)

  mvp_result <- expanding_tvmvp(
    tmp,
    initial_window = 60,
    rebal_period   = 5,
    max_factors    = 10,
    return_type    = "daily",
    kernel_func    = epanechnikov_kernel
  )

  # Check that result is an R6 object of class RollingWindow
  expect_true(inherits(mvp_result, "ExpandingWindow"))

  # Optionally check something else
})
