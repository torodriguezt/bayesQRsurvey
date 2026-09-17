# Tests for posterior_interval(), the rstantools generic used in place of
# confint().

post_fit <- function(quantile = c(0.25, 0.5), ...) {
  set.seed(52)
  n <- 150
  d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), w = runif(n, 1, 3))
  d$y <- 2 + 1.5 * d$x1 - 0.8 * d$x2 + rnorm(n)
  set.seed(53)
  fit <- suppressWarnings(
    bqr.svy(y ~ x1 + x2, weights = w, data = d, quantile = quantile,
            niter = 1000, burnin = 200, verbose = FALSE, ...)
  )
  list(fit = fit, data = d)
}


test_that("posterior_interval returns equal-tailed intervals at the 95% level", {
  fit <- post_fit()$fit

  pi <- posterior_interval(fit, tau = 0.5)
  expect_equal(dim(pi), c(3L, 2L))
  expect_equal(rownames(pi), c("(Intercept)", "x1", "x2"))
  # prob defaults to 0.95, matching every other interval the package reports
  expect_equal(colnames(pi), c("2.5%", "97.5%"))
  expect_true(all(pi[, 1] < pi[, 2]))

  # the point estimate lies inside
  cf <- coef(fit, tau = 0.5)
  expect_true(all(pi[, 1] < cf & cf < pi[, 2]))

  # a narrower prob gives a narrower interval
  pi90 <- posterior_interval(fit, tau = 0.5, prob = 0.9)
  expect_equal(colnames(pi90), c("5%", "95%"))
  expect_true(all(pi90[, 2] - pi90[, 1] < pi[, 2] - pi[, 1]))

  # and the default agrees exactly with what summary() reports
  blk <- summary(fit)$per_tau[["tau=0.500"]]$coef_summary
  expect_equal(unname(pi[, 1]), blk$lower_ci)
  expect_equal(unname(pi[, 2]), blk$upper_ci)

  expect_error(posterior_interval(fit, tau = 0.5, prob = 1.5),
               "must be a single number")
  expect_error(posterior_interval(fit, tau = 0.42), "was not fitted")
})


test_that("posterior_interval is the interval method; confint is left to R", {
  # The package defines no confint method: credible intervals are reported by
  # posterior_interval(), and confint() is left to base R rather than redefined.
  set.seed(52)
  n <- 120
  d <- data.frame(x1 = runif(n, -1, 1), w = runif(n, 1, 3))
  d$y <- 1 + d$x1 + rnorm(n)
  set.seed(53)
  fit <- suppressWarnings(
    bqr.svy(y ~ x1, weights = w, data = d, quantile = 0.5,
            niter = 600, burnin = 200, verbose = FALSE))

  expect_null(getS3method("confint", "bqr.svy", optional = TRUE))

  ci <- posterior_interval(fit, prob = 0.95)
  expect_true(is.matrix(ci))
  expect_equal(rownames(ci), names(coef(fit)))

  # it reports quantiles of the draws, which is what summary() shows
  D <- as.matrix(fit, tau = 0.5)
  expect_equal(unname(ci[, 1]), unname(apply(D, 2, quantile, probs = 0.025)))
})
