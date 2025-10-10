# =====================================================
# Tests for plotting utilities
# =====================================================
test_that("plot.bqr.svy returns a ggplot object", {
  set.seed(131415)
  data <- data.frame(x = seq(-2, 2, length.out = 20),
                     y = 1 + 0.5 * seq(-2, 2, length.out = 20) + rnorm(20, 0, 0.3))

  fit <- bqr.svy(y ~ x, data = data, quantile = 0.5, niter = 500)
  result <- plot(fit, which = "x1", tau = 0.5)

  expect_s3_class(result, "ggplot")
  expect_true(is.data.frame(result$data))
  expect_true(all(c("x", "y") %in% names(result$data)))
})

