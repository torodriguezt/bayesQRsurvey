# Tests for the mo.bqr.svy extractor methods: coef(), sigma() and nobs().
#
# Estimation is by EM, so these are point summaries. There is deliberately no
# vcov() or posterior_interval() for this class; the last block pins that down.

mo_fit <- function(quantile = c(0.25, 0.5), n_dir = 4, ...) {
  set.seed(61)
  n <- 150
  d <- data.frame(x1 = runif(n, -1, 1), w = runif(n, 1, 3))
  d$y1 <- 1 + 1.2 * d$x1 + rnorm(n)
  d$y2 <- 0.5 - 0.7 * d$x1 + rnorm(n)
  set.seed(62)
  fit <- suppressWarnings(
    mo.bqr.svy(cbind(y1, y2) ~ x1, weights = w, data = d, quantile = quantile,
               n_dir = n_dir, max_iter = 500, verbose = FALSE, ...)
  )
  list(fit = fit, data = d)
}


test_that("coef returns parameters by directions", {
  fit <- mo_fit()$fit

  cf <- coef(fit, tau = 0.5)
  expect_true(is.matrix(cf))
  expect_equal(rownames(cf), c("(Intercept)", "x1", "gamma_1"))
  expect_equal(colnames(cf), paste0("dir_", 1:4))
  expect_true(all(is.finite(cf)))

  # one entry per quantile when tau is not given
  cfa <- coef(fit)
  expect_length(cfa, 2L)
  expect_named(cfa, c("tau=0.250", "tau=0.500"))
  expect_equal(cfa[["tau=0.500"]], cf)
})


test_that("direction selects a subset of directions", {
  fit <- mo_fit()$fit

  expect_equal(colnames(coef(fit, tau = 0.5, direction = 1:2)),
               c("dir_1", "dir_2"))
  expect_equal(coef(fit, tau = 0.5, direction = 3),
               coef(fit, tau = 0.5)[, "dir_3", drop = FALSE])
  expect_length(sigma(fit, tau = 0.5, direction = c(1, 4)), 2L)

  expect_error(coef(fit, tau = 0.5, direction = 9), "Invalid direction")
  expect_error(coef(fit, tau = 0.5, direction = 0), "Invalid direction")
  expect_error(coef(fit, tau = 0.5, direction = "dir_1"),
               "must be a numeric vector")
})


test_that("sigma reports the scale for each direction", {
  # held at 1 unless estimated
  fit <- mo_fit()$fit
  s <- sigma(fit, tau = 0.5)
  expect_length(s, 4L)
  expect_named(s, paste0("dir_", 1:4))
  expect_equal(unname(s), rep(1, 4))

  # estimated: a positive value per direction, no longer all ones
  est <- mo_fit(quantile = 0.5, estimate_sigma = TRUE)$fit
  se <- sigma(est)
  expect_length(se, 4L)
  expect_true(all(se > 0))
  expect_false(isTRUE(all.equal(unname(se), rep(1, 4))))

  # one entry per quantile when tau is not given
  expect_length(sigma(fit), 2L)
})


test_that("quantiles are selected numerically, as elsewhere in the package", {
  fit <- mo_fit()$fit
  one <- mo_fit(quantile = 0.5)$fit

  # a single-quantile fit needs no tau
  expect_true(is.matrix(coef(one)))

  expect_error(coef(fit, tau = 0.42), "was not fitted")
  expect_error(sigma(fit, tau = 0.42), "was not fitted")
  expect_error(coef(fit, tau = "0.5"), "must be a numeric vector")
})


test_that("nobs reports the number of observations", {
  obj <- mo_fit()
  expect_identical(nobs(obj$fit), nrow(obj$data))
  expect_type(nobs(obj$fit), "integer")
})


test_that("weights returns the survey weights as supplied, not the normalised copy", {
  obj <- mo_fit()
  w <- weights(obj$fit)
  expect_type(w, "double")
  expect_length(w, nobs(obj$fit))
  expect_equal(w, obj$data$w)
  # the EM works on weights scaled to mean 1; the extractor must not return those
  expect_false(isTRUE(all.equal(mean(w), 1)))
})


test_that("the EM fit offers no posterior uncertainty, and vcov says so", {
  # The object stores posterior modes and no draws, so a covariance matrix
  # cannot be computed. vcov() is defined anyway, so that the generic fails with
  # that reason instead of reporting a missing method.
  obj <- mo_fit()

  expect_true(is.function(getS3method("vcov", "mo.bqr.svy", optional = TRUE)))
  expect_error(vcov(obj$fit), "posterior mode")
  expect_error(vcov(obj$fit), "covariance matrix")

  # neither interval generic is registered for this class: it cannot compute a
  # credible interval, and the package does not redefine confint()
  expect_null(getS3method("posterior_interval", "mo.bqr.svy", optional = TRUE))
  expect_null(getS3method("confint", "mo.bqr.svy", optional = TRUE))
})


test_that("plot draws quantile regions using the data stored in the fit", {
  skip_if_not_installed("ggplot2")
  obj <- mo_fit(quantile = c(0.1, 0.25), n_dir = 8)
  fit <- obj$fit

  # the response names and the data come from the object; this model has one
  # covariate plus the intercept
  r <- plot(fit, xValue = c(1, 0), ngridpoints = 40)
  expect_s3_class(r$plot, "ggplot")
  expect_true(is.data.frame(r$data))
  expect_true(all(c("tau") %in% names(r$data)))

  # overriding the defaults still works
  r2 <- plot(fit, response = c("y1", "y2"), datafile = obj$data,
             xValue = c(1, 0), ngridpoints = 40, paintedArea = FALSE)
  expect_s3_class(r2$plot, "ggplot")

  # a quantile region is conditional on the covariates, so a mismatched xValue
  # is an error rather than an empty region
  expect_error(plot(fit, ngridpoints = 40),
               "must give one value per covariate")
  expect_error(plot(fit, xValue = c(1, 2, 3), ngridpoints = 40),
               "3 supplied")
})


test_that("plotQuantileRegion is gone in favour of the plot method", {
  expect_false("plotQuantileRegion" %in% getNamespaceExports("bayesQRsurvey"))
  expect_false(is.null(getS3method("plot", "mo.bqr.svy", optional = TRUE)))
})


test_that("the model specification methods work for the multiple-output fit", {
  obj <- mo_fit(quantile = 0.5)
  fit <- obj$fit

  expect_equal(formula(fit), cbind(y1, y2) ~ x1, ignore_attr = TRUE)
  expect_s3_class(terms(fit), "terms")

  # the design matrix holds the covariates only; the response is the cbind()
  X <- model.matrix(fit)
  expect_equal(dim(X), c(nrow(obj$data), 2L))
  expect_equal(colnames(X), c("(Intercept)", "x1"))

  # update() evaluates the stored call in the caller, so the data has to be in
  # scope there -- the same requirement lm() has. Build the fit here rather than
  # through mo_fit(), whose data frame is local to that helper.
  d <- obj$data
  local_fit <- suppressWarnings(
    mo.bqr.svy(cbind(y1, y2) ~ x1, weights = w, data = d, quantile = 0.5,
               n_dir = 4, max_iter = 500, verbose = FALSE)
  )
  cl <- update(local_fit, quantile = 0.25, evaluate = FALSE)
  expect_true(is.call(cl))

  refit <- suppressWarnings(update(local_fit, quantile = 0.25))
  expect_s3_class(refit, "mo.bqr.svy")
  expect_equal(refit$quantile, 0.25)
})
