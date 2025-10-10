# Test for multiple quantile functionality
test_that("mo.bqr.svy basic functionality works", {
  set.seed(789)  # Reproducibility

  # Generate test data
  n <- 25  # Small sample for speed
  x <- rnorm(n)
  y <- 1 + 2*x + rnorm(n, 0, 0.5)
  data <- data.frame(x = x, y = y)

  # Test multiple quantile regression
  quantiles <- c(0.25, 0.5, 0.75)
  fit <- mo.bqr.svy(y ~ x, data = data, quantile = quantiles, max_iter = 1500)

  # Test basic structure
  expect_s3_class(fit, "mo.bqr.svy")
  expect_true(is.list(fit))
  expect_true("fit" %in% names(fit))
  expect_true("quantile" %in% names(fit))
  expect_true("call" %in% names(fit))

  # Test dimensions
  expect_equal(length(fit$fit), 3)  # three quantiles
  expect_equal(fit$quantile, quantiles)

  # Test that each fit has proper structure
  for (i in 1:3) {
    expect_true("directions" %in% names(fit$fit[[i]]))
    expect_true("dir_1" %in% names(fit$fit[[i]]$directions))
    expect_true("beta" %in% names(fit$fit[[i]]$directions$dir_1))
    expect_true("sigma" %in% names(fit$fit[[i]]$directions$dir_1))

    expect_true(all(is.finite(fit$fit[[i]]$directions$dir_1$beta)))
    expect_true(is.finite(fit$fit[[i]]$directions$dir_1$sigma))
  }


  # Extract slope coefficients (assuming "x" is the 2nd coefficient)
  beta1_values <- sapply(fit$fit, function(f) {
    f$directions$dir_1$beta[2]
  })

  # Allow some flexibility due to randomness in simulation
  expect_true(all(diff(beta1_values) >= -0.5))

})

test_that("mo.bqr.svy handles different numbers of quantiles", {
  set.seed(101112)

  n <- 20
  x <- rnorm(n)
  y <- x + rnorm(n, 0, 0.3)
  data <- data.frame(x = x, y = y)

  # Test with 2 quantiles
  fit_2 <- mo.bqr.svy(y ~ x, data = data, quantile = c(0.3, 0.7), max_iter = 1500)
  expect_s3_class(fit_2, "mo.bqr.svy")
  expect_equal(length(fit_2$fit), 2)

  # Test with 5 quantiles
  fit_5 <- mo.bqr.svy(y ~ x, data = data,
                      quantile = c(0.1, 0.25, 0.5, 0.75, 0.9), max_iter = 1500)
  expect_s3_class(fit_5, "mo.bqr.svy")
  expect_equal(length(fit_5$fit), 5)
})

test_that("mo.bqr.svy input validation", {
  n <- 15
  data <- data.frame(x = rnorm(n), y = rnorm(n))

  # Test invalid quantiles
  expect_error(mo.bqr.svy(y ~ x, data = data, quantile = c(0.5, 1.5)))
  expect_error(mo.bqr.svy(y ~ x, data = data, quantile = c(-0.1, 0.5)))
  expect_error(mo.bqr.svy(y ~ x, data = data, quantile = c()))  # empty

  # Test unsorted quantiles
  fit <- mo.bqr.svy(y ~ x, data = data, quantile = c(0.7, 0.3, 0.5), max_iter = 1500)
  expect_equal(fit$quantile, c(0.3, 0.5, 0.7))  # Should be sorted
})

