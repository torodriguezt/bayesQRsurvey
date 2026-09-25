prediction_regression_data <- function() {
  set.seed(401)
  n <- 96L
  x <- rnorm(n)
  z <- runif(n, -1, 1)
  data.frame(y = 1 + 0.7 * x - 0.4 * z + rnorm(n),
             x = x, z = z, sigma = x, sigma.1 = z,
             g = factor(rep(c("a", "b", "c"), length.out = n)),
             w = runif(n, 2, 4), o = seq_len(n) / n)
}

test_that("a coefficient named sigma remains distinct from the ALD scale", {
  d <- prediction_regression_data()
  cases <- list(list(method = "ald"),
                list(method = "ald", estimate_sigma = TRUE),
                list(method = "score"), list(method = "approximate"))
  for (args in cases) {
    fit_one <- function(formula) {
      set.seed(402)
      suppressWarnings(do.call(bqr.svy, c(list(
        formula = formula, data = d, weights = d$w,
        quantile = c(0.25, 0.75), niter = 500, burnin = 100,
        verbose = FALSE), args)))
    }
    fit <- fit_one(y ~ sigma + sigma.1)
    reference <- fit_one(y ~ x + z)
    expected_names <- c("(Intercept)", "sigma", "sigma.1")

    for (tau in c(0.25, 0.75)) {
      D <- as.matrix(fit, tau = tau, include_sigma = FALSE)
      ref_D <- as.matrix(reference, tau = tau, include_sigma = FALSE)
      expect_identical(colnames(D), expected_names)
      expect_equal(unname(D), unname(ref_D))
      expect_equal(unname(vcov(fit, tau)), unname(vcov(reference, tau)))
      expect_equal(unname(posterior_interval(fit, tau = tau)),
                   unname(posterior_interval(reference, tau = tau)))
      expect_equal(sigma(fit, tau), sigma(reference, tau))

      nd <- data.frame(sigma = c(-0.5, 0.5), sigma.1 = 0)
      ref_nd <- data.frame(x = nd$sigma, z = nd$sigma.1)
      expect_equal(predict(fit, nd, tau = tau, interval = "credible"),
                   predict(reference, ref_nd, tau = tau, interval = "credible"))

      all_D <- as.matrix(fit, tau = tau, include_sigma = TRUE)
      expect_identical(anyDuplicated(colnames(all_D)), 0L)
      has_scale <- identical(args$method, "ald")
      expect_equal(ncol(all_D), 3L + as.integer(has_scale))
      if (has_scale) expect_identical(colnames(all_D)[4L], "sigma.2")

      diag <- diagnostics(fit, tau)
      ref_diag <- diagnostics(reference, tau)
      expect_identical(diag$variable[1:3], expected_names)
      expect_equal(diag[, -1], ref_diag[, -1])
      expect_identical(anyDuplicated(diag$variable), 0L)
    }

    summ <- summary(fit)$per_tau[[1L]]$coef_summary
    ref_summ <- summary(reference)$per_tau[[1L]]$coef_summary
    expect_identical(summ$variable[1:3], expected_names)
    expect_equal(summ[, -1], ref_summ[, -1])
    expect_identical(anyDuplicated(summ$variable), 0L)
    expect_output(print(fit), "sigma")

    # Plot data must also retain the coefficient rather than use the scale.
    tr <- plot(fit, type = "trace", which = "sigma", tau = 0.25)
    expect_equal(tr$data$value,
                 as.matrix(reference, tau = 0.25, include_sigma = FALSE)[, "x"])
    den <- plot(fit, type = "density", which = "sigma", tau = 0.25)
    expect_equal(den$data$x, tr$data$value)
    qp <- plot(fit, type = "quantile", which = "sigma", add_ols = FALSE)
    expect_equal(qp$data$est, as.numeric(coef(reference)["x", ]))
  }
})

test_that("stored contrasts govern fitted values and plotting after options change", {
  d <- prediction_regression_data()
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  cases <- list(list(method = "ald"),
                list(method = "ald", estimate_sigma = TRUE),
                list(method = "score"), list(method = "approximate"))
  for (args in cases) {
    options(contrasts = c("contr.sum", "contr.poly"))
    set.seed(403)
    fit <- suppressWarnings(do.call(bqr.svy, c(list(
      formula = y ~ g + x, data = d, weights = d$w,
      quantile = 0.5, niter = 500, burnin = 100, verbose = FALSE), args)))
    expected <- predict(fit, interval = "credible")
    expected_X <- model.matrix(fit)
    curve <- plot(fit, type = "fit", which = "x", at = list(g = "b"))

    options(contrasts = c("contr.helmert", "contr.poly"))
    expect_equal(model.matrix(fit), expected_X)
    expect_equal(fitted(fit), expected[, "fit"])
    expect_equal(predict(fit, interval = "credible"), expected)
    expect_equal(predict(fit, newdata = d, interval = "credible"), expected)
    expect_equal(plot(fit, type = "fit", which = "x", at = list(g = "b"))$data,
                 curve$data)
  }
})

test_that("unsupported offsets and invalid weights fail before sampling", {
  d <- prediction_regression_data()
  for (method in c("ald", "score", "approximate")) {
    expect_error(bqr.svy(y ~ x + offset(o), data = d, method = method,
                         niter = 20, verbose = FALSE), "offset.*not supported")
    for (bad in c(0, -1, NA_real_, NaN, Inf, -Inf)) {
      invalid_weights <- d$w
      invalid_weights[1L] <- bad
      expect_error(bqr.svy(y ~ x, data = d, weights = invalid_weights, method = method,
                           niter = 20, verbose = FALSE), "finite.*positive")
    }
  }
})
