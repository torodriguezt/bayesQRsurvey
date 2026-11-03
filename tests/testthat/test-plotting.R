# =====================================================
# Tests for plotting utilities
# =====================================================
test_that("plot.bqr.svy returns a ggplot object", {
  set.seed(131415)
  data <- data.frame(x = seq(-2, 2, length.out = 20),
                     y = 1 + 0.5 * seq(-2, 2, length.out = 20) + rnorm(20, 0, 0.3))

  fit <- bqr.svy(y ~ x, data = data, quantile = 0.5, niter = 500)
  result <- plot(fit, type = "fit", which = "x", tau = 0.5)

  expect_s3_class(result, "ggplot")
  expect_true(is.data.frame(result$data))
  expect_true(all(c("x", "y") %in% names(result$data)))
})


test_that("the fitted curve drawn by plot() is x' beta_bar over the same grid", {
  set.seed(2718)
  n <- 120
  d <- data.frame(x1 = runif(n, -1, 1), g = factor(sample(c("a", "b"), n, TRUE)),
                  w = runif(n, 1, 3))
  d$y <- 1 + 2 * d$x1 - 0.5 * d$x1^2 + (d$g == "b") + rnorm(n)

  set.seed(2719)
  fit <- suppressWarnings(
    bqr.svy(y ~ x1 + I(x1^2) + g, weights = w, data = d,
            quantile = c(0.25, 0.75), niter = 1000, burnin = 200, verbose = FALSE)
  )

  p <- plot(fit, type = "fit", which = "x1", grid_length = 25,
            at = list(g = "b"), show_ci = TRUE)
  pd <- p$data

  for (tt in c(0.25, 0.75)) {
    rows <- pd[pd$tau_numeric == tt, ]
    # the design matrix of the grid at g = "b": (Intercept), x1, I(x1^2), gb
    X    <- cbind(1, rows$x, rows$x^2, 1)
    XB   <- X %*% t(as.matrix(fit, tau = tt))
    band <- t(apply(XB, 1, quantile, probs = c(0.025, 0.975)))
    # the curve is the conditional quantile at the posterior mean, x' beta_bar
    expect_equal(unname(rows$y), unname(drop(X %*% coef(fit, tau = tt))))
    # and the band is the row quantiles of x' beta over the draws
    expect_equal(unname(rows$y_lower), unname(band[, 1]))
    expect_equal(unname(rows$y_upper), unname(band[, 2]))
  }
})

