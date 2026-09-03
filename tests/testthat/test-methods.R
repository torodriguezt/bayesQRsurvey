# Tests for coef(), fitted() and confint().
#
# The point is to notice if something breaks: shapes, names, quantile selection
# and agreement with summary(). Nothing here depends on exact numeric output.

fit_small <- function(quantile = c(0.25, 0.5, 0.75), ...) {
  set.seed(42)
  n <- 120
  d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), w = runif(n, 1, 3))
  d$y <- 2 + 1.5 * d$x1 - 0.8 * d$x2 + rnorm(n)
  set.seed(43)
  fit <- suppressWarnings(
    bqr.svy(y ~ x1 + x2, weights = w, data = d, quantile = quantile,
            niter = 1000, burnin = 200, verbose = FALSE, ...)
  )
  list(fit = fit, data = d)
}


test_that("coef gives a vector for one quantile and a matrix for several", {
  one   <- fit_small(quantile = 0.5)$fit
  three <- fit_small()$fit

  cf1 <- coef(one)
  expect_null(dim(cf1))
  expect_named(cf1, c("(Intercept)", "x1", "x2"))
  expect_true(all(is.finite(cf1)))

  cf3 <- coef(three)
  expect_equal(dim(cf3), c(3L, 3L))
  expect_equal(rownames(cf3), c("(Intercept)", "x1", "x2"))
  expect_equal(colnames(cf3), c("tau=0.250", "tau=0.500", "tau=0.750"))
  expect_true(all(is.finite(cf3)))

  # picking one quantile out of several gives the matching column
  expect_equal(coef(three, tau = 0.5), cf3[, "tau=0.500"])
})


test_that("fitted has one value per observation", {
  obj <- fit_small()
  fit <- obj$fit

  f1 <- fitted(fit, tau = 0.5)
  expect_length(f1, nrow(obj$data))
  expect_true(all(is.finite(f1)))

  # fitted is X %*% coef
  X <- model.matrix(y ~ x1 + x2, data = obj$data)
  expect_equal(unname(f1), as.numeric(X %*% coef(fit, tau = 0.5)))

  # one column per quantile when tau is not given
  expect_equal(dim(fitted(fit)), c(nrow(obj$data), 3L))
})


test_that("confint returns sensible intervals that agree with summary", {
  fit <- fit_small()$fit

  ci <- confint(fit, tau = 0.5)
  expect_equal(dim(ci), c(3L, 2L))
  expect_equal(rownames(ci), c("(Intercept)", "x1", "x2"))
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))

  # lower < upper, and the point estimate sits inside
  expect_true(all(ci[, 1] < ci[, 2]))
  cf <- coef(fit, tau = 0.5)
  expect_true(all(ci[, 1] < cf & cf < ci[, 2]))

  # identical to what summary() reports, by construction
  blk <- summary(fit)$per_tau[["tau=0.500"]]$coef_summary
  expect_equal(unname(ci[, 1]), blk$lower_ci)
  expect_equal(unname(ci[, 2]), blk$upper_ci)

  # a lower level gives a narrower interval
  ci80 <- confint(fit, tau = 0.5, level = 0.80)
  expect_true(all(ci80[, 2] - ci80[, 1] < ci[, 2] - ci[, 1]))

  # one entry per quantile when tau is not given
  expect_length(confint(fit), 3L)
})


test_that("parm selects coefficients by name or position", {
  fit <- fit_small()$fit

  expect_equal(rownames(confint(fit, parm = "x1", tau = 0.5)), "x1")
  expect_equal(rownames(confint(fit, parm = 2, tau = 0.5)), "x1")
  expect_equal(rownames(confint(fit, parm = c("x1", "x2"), tau = 0.5)),
               c("x1", "x2"))
})


test_that("quantiles are selected numerically, with clear errors", {
  fit  <- fit_small()$fit
  one  <- fit_small(quantile = 0.5)$fit

  # a single-quantile fit needs no tau at all
  expect_length(coef(one), 3L)

  expect_error(coef(fit, tau = 0.42), "was not fitted")
  expect_error(fitted(fit, tau = 0.42), "was not fitted")
  expect_error(confint(fit, tau = 0.42), "was not fitted")
  expect_error(coef(fit, tau = "0.5"), "must be a numeric vector")
  expect_error(confint(fit, tau = 0.5, level = 1.5), "must be a single number")
  expect_error(confint(fit, parm = "nope", tau = 0.5), "Unknown parameter")
})


test_that("the methods work for all three estimation methods", {
  for (m in c("ald", "score", "approximate")) {
    fit <- fit_small(quantile = c(0.25, 0.75), method = m)$fit

    expect_true(all(is.finite(coef(fit))), info = m)
    expect_true(all(is.finite(fitted(fit))), info = m)
    expect_true(all(is.finite(confint(fit, tau = 0.25))), info = m)
  }
})


test_that("diagnostics replaces reaching into $diagnosis", {
  fit <- fit_small()$fit

  d <- diagnostics(fit, tau = 0.5)
  expect_s3_class(d, "bqr.svy.diagnostics")
  expect_equal(names(d), c("variable", "rhat", "ess_bulk", "ess_tail"))
  expect_equal(d$variable, c("(Intercept)", "x1", "x2"))

  # same content the internal element holds, without the label gymnastics
  expect_equal(d, fit$diagnosis[["tau=0.500"]])

  expect_length(diagnostics(fit), 3L)
  expect_error(diagnostics(fit, tau = 0.42), "was not fitted")
})
