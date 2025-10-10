# Test for basic bqr.svy functionality
test_that("bqr.svy basic functionality works", {
  set.seed(123)  # Reproducibility

  # Generate simple test data
  n <- 30  # Small sample for speed
  x <- rnorm(n)
  y <- 2 + 3*x + rnorm(n, 0, 0.5)
  data <- data.frame(x = x, y = y)
  w <- rep(1, n)  # equal weights

  # Test basic single quantile regression
  fit <- bqr.svy(y ~ x, data = data, quantile = 0.5, niter = 2000)  # Few iterations for speed

  # Test basic structure
  expect_s3_class(fit, "bqr.svy")
  expect_true(is.list(fit))
  expect_true("beta" %in% names(fit))
  expect_true("quantile" %in% names(fit))
  expect_true("call" %in% names(fit))

  # Test dimensions
  expect_equal(length(fit$beta), 2)  # intercept + x
  expect_equal(fit$quantile, 0.5)

  # Test that coefficients are reasonable (not NaN/Inf)
  expect_true(all(is.finite(fit$beta)))

  # Test different quantiles
  fit_25 <- bqr.svy(y ~ x, data = data, quantile = 0.25, niter = 2000)
  expect_equal(fit_25$quantile, 0.25)
  expect_true(all(is.finite(fit_25$beta)))

  fit_75 <- bqr.svy(y ~ x, data = data, quantile = 0.75, niter = 2000)
  expect_equal(fit_75$quantile, 0.75)
  expect_true(all(is.finite(fit_75$beta)))
})

test_that("bqr.svy handles edge cases", {
  set.seed(456)

  # Test with minimal data
  n <- 10
  x <- 1:n
  y <- x + rnorm(n, 0, 0.1)
  data <- data.frame(x = x, y = y)

  fit <- bqr.svy(y ~ x, data = data, quantile = 0.5, niter = 2000)
  expect_s3_class(fit, "bqr.svy")
  expect_true(all(is.finite(fit$beta)))

  # Test with weights
  w <- runif(n, 0.5, 2)
  fit_w <- bqr.svy(y ~ x, data = data, weights = w, quantile = 0.5, niter = 2000)
  expect_s3_class(fit_w, "bqr.svy")
  expect_true(all(is.finite(fit_w$beta)))
})

test_that("bqr.svy input validation", {
  # Test invalid quantile
  n <- 20
  data <- data.frame(x = rnorm(n), y = rnorm(n))

  expect_error(bqr.svy(y ~ x, data = data, quantile = 0))
  expect_error(bqr.svy(y ~ x, data = data, quantile = 1))
  expect_error(bqr.svy(y ~ x, data = data, quantile = 1.5))
  expect_error(bqr.svy(y ~ x, data = data, quantile = -0.1))

  # Test invalid iterations
  expect_error(bqr.svy(y ~ x, data = data, niter = 0))
  expect_error(bqr.svy(y ~ x, data = data, niter = -5))

  # Test missing data handling
  data_na <- data.frame(x = c(1, 2, NA), y = c(1, 2, 3))
  expect_error(bqr.svy(y ~ x, data = data_na, quantile = 0.5))
})

