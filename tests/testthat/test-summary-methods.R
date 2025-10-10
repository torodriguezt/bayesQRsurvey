# =====================================================
# Tests for summary.bayesQRsurvey (bqr.svy and mo.bqr.svy)
# =====================================================

test_that("summary.bqr.svy returns valid output with realistic sample weights", {
  set.seed(123)

  # --- Generate finite population and sample ---
  N <- 10000
  x1_p <- runif(N, -1, 1)
  x2_p <- runif(N, -1, 1)
  y_p  <- 2 + 1.5 * x1_p - 0.8 * x2_p + rnorm(N)
  n <- 300
  z_aux <- rnorm(N, mean = 1 + y_p, sd = 0.5)
  p_aux <- 1 / (1 + exp(2.5 - 0.5 * z_aux))
  s_ind <- sample(1:N, n, replace = FALSE, prob = p_aux)
  data <- data.frame(
    y  = y_p[s_ind],
    x1 = x1_p[s_ind],
    x2 = x2_p[s_ind],
    w  = 1 / p_aux[s_ind]
  )

  # --- Fit model ---
  fit <- bqr.svy(y ~ x1 + x2, weights = ~w, data = data, quantile = 0.5, niter = 300)
  s <- summary(fit)

  # --- Structure checks ---
  expect_s3_class(s, "summary.bqr.svy")
  expect_true(is.list(s))
  expect_true(length(s) >= 1)
})


# -----------------------------------------------------

test_that("summary.mo.bqr.svy returns valid output with informative prior", {
  set.seed(456)

  # --- Generate population and sample ---
  N <- 5000
  x1_p <- runif(N, -1, 1)
  x2_p <- runif(N, -1, 1)
  y1_p <- 2 + 1.5 * x1_p - 0.8 * x2_p + rnorm(N)
  y2_p <- 1 + 0.5 * x1_p - 0.2 * x2_p + rnorm(N)
  n <- 200
  z_aux <- rnorm(N, mean = 1 + y1_p, sd = 0.5)
  p_aux <- 1 / (1 + exp(2.5 - 0.5 * z_aux))
  s_ind <- sample(1:N, n, replace = FALSE, prob = p_aux)
  data <- data.frame(
    y1 = y1_p[s_ind],
    y2 = y2_p[s_ind],
    x1 = x1_p[s_ind],
    x2 = x2_p[s_ind],
    w  = 1 / p_aux[s_ind]
  )

  # --- Define prior ---
  prior_general <- prior(
    beta_x_mean = c(2, 1.5, -0.8),
    beta_x_cov  = diag(c(0.25, 0.25, 0.25)),
    sigma_shape = 3,
    sigma_rate  = 2,
    beta_y_mean = 1,
    beta_y_cov  = 0.25
  )

  # --- Fit model ---
  fit_mo <- mo.bqr.svy(cbind(y1, y2) ~ x1 + x2,
                       weights = data$w,
                       data = data,
                       prior = prior_general,
                       n_dir = 5, estimate_sigma = TRUE)
  s_mo <- summary(fit_mo)

  # --- Structure checks ---
  expect_s3_class(s_mo, "summary.mo.bqr.svy")
  expect_true(is.list(s_mo))
})
