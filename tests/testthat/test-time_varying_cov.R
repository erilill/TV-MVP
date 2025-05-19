test_that("time_varying_cov returns an appropriate R6 object with covariance estimates", {
  set.seed(123)
  T <- 100  # Number of time periods
  p <- 20   # Number of assets
  returns <- matrix(rnorm(T * p, mean = 0.001, sd = 0.02), ncol = p)

  tmp <- TVMVP$new()
  tmp$set(data = returns)

  cov_mat <- time_varying_cov(tmp,
                              m=1,
                              kernel_func = epanechnikov_kernel,
                              M0 = 10,
                              rho_grid = seq (0.005 , 2,
                                              length.out = 30),
                              floor_value = 1e-12,
                              epsilon2 = 1e-6,
                              full_output = FALSE)

  # R6 class check
  expect_true(is.matrix(cov_mat))

  # Optional: check for something else
})
