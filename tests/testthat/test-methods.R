# Tests for the extractor methods: coef(), fitted(), vcov() and nobs().
#
# The point is to notice if something breaks: shapes, names, quantile selection
# and internal consistency between the methods. Nothing here depends on exact
# numeric output.

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


test_that("quantiles are selected numerically, with clear errors", {
  fit  <- fit_small()$fit
  one  <- fit_small(quantile = 0.5)$fit

  # a single-quantile fit needs no tau at all
  expect_length(coef(one), 3L)

  expect_error(coef(fit, tau = 0.42), "was not fitted")
  expect_error(fitted(fit, tau = 0.42), "was not fitted")
  expect_error(coef(fit, tau = "0.5"), "must be a numeric vector")
})


test_that("the methods work for all three estimation methods", {
  for (m in c("ald", "score", "approximate")) {
    fit <- fit_small(quantile = c(0.25, 0.75), method = m)$fit

    expect_true(all(is.finite(coef(fit))), info = m)
    expect_true(all(is.finite(fitted(fit))), info = m)
    expect_true(all(is.finite(vcov(fit, tau = 0.25))), info = m)
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


test_that("vcov is the posterior covariance of the coefficient draws", {
  fit <- fit_small()$fit

  V <- vcov(fit, tau = 0.5)
  expect_equal(dim(V), c(3L, 3L))
  expect_equal(rownames(V), c("(Intercept)", "x1", "x2"))
  expect_true(isSymmetric(V))
  expect_true(all(diag(V) > 0))

  # it is exactly var() of the posterior draws; the class only governs printing
  D <- bayesQRsurvey:::.draws_for_tau(fit, 0.5)
  expect_equal(unclass(V), stats::var(D))
  expect_s3_class(V, "bqr.svy.vcov")
  expect_true(is.matrix(V))

  # the implied standard errors match the spread of the draws
  expect_equal(sqrt(diag(V)), apply(D, 2, stats::sd))

  expect_length(vcov(fit), 3L)
})


test_that("nobs reports the number of observations", {
  obj <- fit_small()
  expect_identical(nobs(obj$fit), nrow(obj$data))
  expect_type(nobs(obj$fit), "integer")
})


test_that("weights returns the survey weights as supplied", {
  obj <- fit_small()
  w <- weights(obj$fit)
  expect_type(w, "double")
  expect_length(w, nobs(obj$fit))
  # the weights as given, not a normalised copy
  expect_equal(w, obj$data$w)

  # a fit without weights records unit weights
  unw <- suppressWarnings(
    bqr.svy(y ~ x1 + x2, data = obj$data, quantile = 0.5,
            niter = 500, burnin = 100, verbose = FALSE)
  )
  expect_equal(weights(unw), rep(1, nrow(obj$data)))
})


test_that("as.matrix returns the posterior draws", {
  fit  <- fit_small()$fit
  one  <- fit_small(quantile = 0.5)$fit

  D <- as.matrix(fit, tau = 0.5)
  expect_true(is.matrix(D))
  expect_equal(colnames(D), c("(Intercept)", "x1", "x2"))
  expect_equal(nrow(D), 800L)          # (niter - burnin) / thin
  expect_true(all(is.finite(D)))

  # it is the same content the internal component holds, without the label
  expect_equal(D, fit$draws[["tau=0.500"]][, colnames(D), drop = FALSE])

  # collapsing the draws gives the point summaries
  expect_equal(colMeans(D), coef(fit, tau = 0.5))
  expect_equal(stats::var(D), unclass(vcov(fit, tau = 0.5)))

  # a single-quantile fit needs no tau
  expect_equal(dim(as.matrix(one)), c(800L, 3L))
  expect_error(as.matrix(fit), "specify 'tau'")
  expect_error(as.matrix(fit, tau = 0.42), "was not fitted")
})


test_that("as.matrix supports the inference the package does not precompute", {
  fit <- fit_small()$fit

  # posterior probability
  p <- mean(as.matrix(fit, tau = 0.5)[, "x1"] > 0)
  expect_gte(p, 0)
  expect_lte(p, 1)

  # contrast between two quantiles, the question specific to quantile regression
  d75 <- as.matrix(fit, tau = 0.75)[, "x1"]
  d25 <- as.matrix(fit, tau = 0.25)[, "x1"]
  expect_length(d75, length(d25))
  expect_true(is.finite(mean(d75 - d25)))
})


test_that("as.matrix follows the sigma column convention", {
  plain <- fit_small(quantile = 0.5)$fit
  est   <- fit_small(quantile = 0.5, estimate_sigma = TRUE)$fit

  # the constant sigma column is hidden unless estimated, or asked for
  expect_false("sigma" %in% colnames(as.matrix(plain)))
  expect_true("sigma" %in% colnames(as.matrix(plain, include_sigma = TRUE)))
  expect_true("sigma" %in% colnames(as.matrix(est)))
})


test_that("sigma reports the scale of the working likelihood", {
  # "ald" estimating sigma: the posterior mean of the sigma draws
  est <- fit_small(quantile = 0.5, estimate_sigma = TRUE)$fit
  expect_true(is.finite(sigma(est)))
  expect_equal(sigma(est), mean(as.matrix(est)[, "sigma"]))

  # "ald" holding it fixed: exactly 1, which is what the model uses
  plain <- fit_small(quantile = 0.5)$fit
  expect_identical(sigma(plain), 1)

  # score and approximate have no scale parameter
  for (m in c("score", "approximate")) {
    f <- fit_small(quantile = 0.5, method = m)$fit
    expect_true(is.na(sigma(f)), info = m)
  }

  # one entry per quantile, named, when several were fitted
  many <- fit_small()$fit
  s <- sigma(many)
  expect_length(s, 3L)
  expect_named(s, c("tau=0.250", "tau=0.500", "tau=0.750"))
})


test_that("the model specification methods recover what the fit was built on", {
  obj <- fit_small(quantile = 0.5)
  fit <- obj$fit

  expect_equal(formula(fit), y ~ x1 + x2, ignore_attr = TRUE)
  expect_s3_class(terms(fit), "terms")

  X <- model.matrix(fit)
  expect_equal(dim(X), c(nrow(obj$data), 3L))
  expect_equal(colnames(X), c("(Intercept)", "x1", "x2"))
  # the design matrix must agree with the one fitted() works from
  expect_equal(unname(drop(X %*% coef(fit))), unname(fitted(fit)))
})


test_that("update refits with part of the call changed", {
  d <- fit_small(quantile = 0.5)$data
  fit <- suppressWarnings(bqr.svy(y ~ x1 + x2, weights = w, data = d,
                                  quantile = 0.5, niter = 1000, burnin = 200,
                                  verbose = FALSE))

  # evaluate = FALSE returns the modified call rather than running it
  cl <- update(fit, quantile = 0.25, evaluate = FALSE)
  expect_true(is.call(cl))
  expect_equal(eval(cl$quantile), 0.25)

  # a formula change drops the term and the refit has one coefficient fewer
  smaller <- suppressWarnings(update(fit, . ~ . - x2))
  expect_s3_class(smaller, "bqr.svy")
  expect_named(coef(smaller), c("(Intercept)", "x1"))
  expect_equal(nobs(smaller), nobs(fit))
})


test_that("predict evaluates the conditional quantile at new covariates", {
  obj <- fit_small(quantile = c(0.25, 0.5, 0.75))
  fit <- obj$fit
  nd  <- data.frame(x1 = c(-0.5, 0, 0.5), x2 = 0)

  # without newdata, predict is fitted(), as in predict.rq
  expect_equal(predict(fit, tau = 0.5), fitted(fit, tau = 0.5))

  p1 <- predict(fit, nd, tau = 0.5)
  expect_null(dim(p1))
  expect_length(p1, 3L)
  expect_equal(unname(p1),
               unname(drop(cbind(1, nd$x1, nd$x2) %*% coef(fit, tau = 0.5))))

  # several quantiles give one column each
  p3 <- predict(fit, nd)
  expect_equal(dim(p3), c(3L, 3L))
  expect_equal(colnames(p3), names(diagnostics(fit)))

  # a quantile that was not fitted is an error, never the nearest one
  expect_error(predict(fit, nd, tau = 0.4), "not fitted")
})


test_that("the credible band brackets the fit and agrees with the draws", {
  fit <- fit_small(quantile = 0.5)$fit
  nd  <- data.frame(x1 = c(-0.5, 0.5), x2 = 0)

  ci <- predict(fit, nd, tau = 0.5, interval = "credible")
  expect_equal(colnames(ci), c("fit", "lower", "higher"))
  expect_true(all(ci[, "lower"] <= ci[, "fit"]))
  expect_true(all(ci[, "fit"]   <= ci[, "higher"]))

  # a wider level gives a wider band
  wide <- predict(fit, nd, tau = 0.5, interval = "credible", level = 0.99)
  expect_true(all(wide[, "higher"] - wide[, "lower"] >=
                  ci[, "higher"] - ci[, "lower"]))

  # the band is exactly the row quantiles of X %*% t(draws)
  X  <- cbind(1, nd$x1, nd$x2)
  XB <- X %*% t(as.matrix(fit, tau = 0.5))
  expect_equal(unname(ci[, "lower"]),
               unname(apply(XB, 1, quantile, probs = 0.025)))

  expect_error(predict(fit, nd, interval = "prediction"), "arg")
  expect_error(predict(fit, nd, interval = "credible", level = 1), "level")
})


test_that("predict rebuilds factors and transformed terms from newdata", {
  set.seed(91)
  n <- 200
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)),
                  x = runif(n, 0.1, 5), w = runif(n, 1, 2))
  d$y <- 1 + as.numeric(d$g) + 2 * log(d$x) + rnorm(n)
  fit <- suppressWarnings(
    bqr.svy(y ~ g + log(x), weights = w, data = d, quantile = 0.5,
            niter = 800, burnin = 200, verbose = FALSE)
  )

  # a level absent from newdata must still get the coding of the fit
  nd <- data.frame(g = factor("b", levels = levels(d$g)), x = 2)
  p  <- predict(fit, nd, tau = 0.5)
  expect_length(p, 1L)
  expect_true(is.finite(p))

  cf <- coef(fit, tau = 0.5)
  expect_equal(unname(p),
               unname(cf[["(Intercept)"]] + cf[["gb"]] + cf[["log(x)"]] * log(2)))

  # an unseen factor level is an error, not a silent NA
  bad <- data.frame(g = factor("z"), x = 2)
  expect_error(predict(fit, bad, tau = 0.5))
})


test_that("predict accepts character columns for variables fitted as factors", {
  set.seed(92)
  n <- 200
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)),
                  x = runif(n, -1, 1), w = runif(n, 1, 2))
  d$y <- 1 + as.numeric(d$g) + 2 * d$x + rnorm(n)
  fit <- suppressWarnings(
    bqr.svy(y ~ g + x, weights = w, data = d, quantile = 0.5,
            niter = 800, burnin = 200, verbose = FALSE)
  )

  # data.frame() no longer builds factors from strings, so newdata written the
  # obvious way arrives as character; the stored xlevels must absorb that
  chr <- data.frame(g = c("a", "b"), x = 0)
  fac <- data.frame(g = factor(c("a", "b"), levels = levels(d$g)), x = 0)
  expect_equal(predict(fit, chr, tau = 0.5), predict(fit, fac, tau = 0.5))

  # a single level in newdata keeps the coding of the fit instead of becoming
  # its own baseline
  one <- predict(fit, data.frame(g = "b", x = 0), tau = 0.5)
  cf  <- coef(fit, tau = 0.5)
  expect_equal(unname(one), unname(cf[["(Intercept)"]] + cf[["gb"]]))

  ci <- predict(fit, chr, tau = 0.5, interval = "credible")
  expect_equal(colnames(ci), c("fit", "lower", "higher"))
  expect_true(all(is.finite(ci)))
})
