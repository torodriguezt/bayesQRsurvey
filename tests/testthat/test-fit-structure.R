# Structural tests for the object returned by bqr.svy().
#
# A safety net for the 0.4.0 refactor: the class collapsed to a single one, and
# the object now has the same shape whether one or several quantiles were fitted.
# These are the invariants that must not silently drift.

toy_data <- function(n = 100) {
  set.seed(7)
  d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), w = runif(n, 1, 3))
  d$y <- 2 + 1.5 * d$x1 - 0.8 * d$x2 + rnorm(n)
  d
}

toy_fit <- function(quantile = 0.5, ...) {
  d <- toy_data()
  set.seed(8)
  suppressWarnings(
    bqr.svy(y ~ x1 + x2, weights = w, data = d, quantile = quantile,
            niter = 1000, burnin = 200, verbose = FALSE, ...)
  )
}


test_that("the fitted object has a single class", {
  expect_identical(class(toy_fit(0.5)), "bqr.svy")
  expect_identical(class(toy_fit(c(0.25, 0.75))), "bqr.svy")
})


test_that("the object shape does not depend on the number of quantiles", {
  one   <- toy_fit(0.5)
  three <- toy_fit(c(0.25, 0.5, 0.75))

  for (fit in list(one, three)) {
    expect_true(is.matrix(fit$beta))
    expect_true(is.list(fit$draws))
    expect_true(is.list(fit$diagnosis))
    expect_equal(rownames(fit$beta), c("(Intercept)", "x1", "x2"))
  }

  expect_equal(ncol(one$beta), 1L)
  expect_equal(ncol(three$beta), 3L)

  expect_named(one$draws, "tau=0.500")
  expect_named(three$draws, c("tau=0.250", "tau=0.500", "tau=0.750"))
  expect_named(three$diagnosis, c("tau=0.250", "tau=0.500", "tau=0.750"))
  expect_length(three$accept_rate, 3L)
})


test_that("quantiles are sorted and de-duplicated", {
  fit <- toy_fit(c(0.75, 0.25))
  expect_equal(fit$quantile, c(0.25, 0.75))

  expect_warning(
    dup <- suppressWarnings(toy_fit(c(0.5, 0.5))),
    NA  # toy_fit already suppresses; just check the result
  )
  expect_equal(dup$quantile, 0.5)
})


test_that("niter, burnin and thin control the number of retained draws", {
  d <- toy_data()
  run <- function(method, niter, burnin, thin) {
    set.seed(8)
    f <- suppressWarnings(
      bqr.svy(y ~ x1 + x2, weights = w, data = d, method = method, quantile = 0.5,
              niter = niter, burnin = burnin, thin = thin, verbose = FALSE))
    nrow(f$draws[[1]])
  }

  for (m in c("ald", "score", "approximate")) {
    expect_equal(run(m, 1000,   0, 1), 1000, info = m)
    expect_equal(run(m, 1000, 200, 1),  800, info = m)
    expect_equal(run(m, 1000, 200, 2),  400, info = m)
    expect_equal(run(m, 1000,   0, 5),  200, info = m)
  }
})


test_that("the fit records what it was given", {
  d <- toy_data()
  fit <- toy_fit(c(0.25, 0.75))

  expect_equal(fit$n, nrow(d))
  expect_equal(fit$weights, d$w)
  expect_equal(fit$warmup, 200)
  expect_equal(fit$thin, 1)
  expect_equal(fit$method, "ald")
  expect_equal(deparse(fit$formula), "y ~ x1 + x2")

  # a real call object, so getCall() and update() work
  expect_true(is.call(fit$call))
  expect_true(is.call(getCall(fit)))
})


test_that("the sigma column follows estimate_sigma", {
  plain <- toy_fit(0.5)
  est   <- toy_fit(0.5, estimate_sigma = TRUE)

  # "ald" always returns a sigma column; it is constant when not estimated
  expect_true("sigma" %in% colnames(plain$draws[[1]]))
  expect_equal(unique(plain$draws[[1]][, "sigma"]), 1)
  expect_gt(stats::sd(est$draws[[1]][, "sigma"]), 0)

  # and it is reported by summary() only when it was estimated
  expect_false("sigma" %in% summary(plain)$per_tau[[1]]$coef_summary$variable)
  expect_true("sigma" %in% summary(est)$per_tau[[1]]$coef_summary$variable)
})


test_that("the package defines no method on the base list class", {
  expect_false("summary.list" %in% ls(asNamespace("bayesQRsurvey")))
  expect_s3_class(summary(list(a = 1, b = "x")), "table")
})
