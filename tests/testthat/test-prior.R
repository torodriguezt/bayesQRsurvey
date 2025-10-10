# =====================================================
# Tests for prior() constructor
# =====================================================

test_that("prior() builds a valid prior object with all parameters", {
  beta_x_mean <- c(2, 1.5, -0.8)
  beta_x_cov  <- diag(c(0.25, 0.25, 0.25))
  beta_y_mean <- 1
  beta_y_cov  <- 0.25

  prior_general <- prior(
    beta_x_mean = beta_x_mean,
    beta_x_cov  = beta_x_cov,
    sigma_shape = 3,
    sigma_rate  = 2,
    beta_y_mean = beta_y_mean,
    beta_y_cov  = beta_y_cov
  )

  # --- Structure and class ---
  expect_s3_class(prior_general, "prior")
  expect_true(is.list(prior_general))

  # --- Core names ---
  expect_true(all(c("beta_x_mean", "beta_x_cov",
                    "sigma_shape", "sigma_rate",
                    "beta_y_mean", "beta_y_cov") %in% names(prior_general)))

  # --- Type checks ---
  expect_type(prior_general$beta_x_mean, "double")
  expect_true(is.matrix(prior_general$beta_x_cov))
  expect_true(is.numeric(prior_general$sigma_shape))
  expect_true(is.numeric(prior_general$sigma_rate))
  expect_type(prior_general$beta_y_mean, "double")
  expect_true(is.numeric(prior_general$beta_y_cov))

  # --- Dimension consistency ---
  expect_equal(length(prior_general$beta_x_mean), nrow(prior_general$beta_x_cov))
  expect_equal(ncol(prior_general$beta_x_cov), length(prior_general$beta_x_mean))
})

# -----------------------------------------------------

test_that("prior() works without optional sigma_* and beta_y_* parameters", {
  beta_x_mean <- c(0.5, -0.2)
  beta_x_cov  <- diag(2)

  prior_simple <- prior(
    beta_x_mean = beta_x_mean,
    beta_x_cov  = beta_x_cov
  )

  expect_s3_class(prior_simple, "prior")
  expect_true(is.list(prior_simple))

  # Should still include sigma defaults even if not provided
  expect_equal(prior_simple$sigma_shape, 0.001)
  expect_equal(prior_simple$sigma_rate, 0.001)

  # Optional elements (beta_y_mean / beta_y_cov) should exist or be NULL
  expect_true("beta_y_mean" %in% names(prior_simple))
  expect_true("beta_y_cov" %in% names(prior_simple))
  expect_true(is.null(prior_simple$beta_y_mean) || is.numeric(prior_simple$beta_y_mean))
  expect_true(is.null(prior_simple$beta_y_cov) || is.numeric(prior_simple$beta_y_cov))
})

# -----------------------------------------------------

test_that("prior() returns sensible defaults when called with no arguments", {
  p_default <- prior()

  expect_s3_class(p_default, "prior")
  expect_true(is.list(p_default))

  # Default hyperparameters should be set
  expect_equal(p_default$sigma_shape, 0.001)
  expect_equal(p_default$sigma_rate, 0.001)

  # Optional arguments should be NULL
  expect_null(p_default$beta_x_mean)
  expect_null(p_default$beta_x_cov)
  expect_null(p_default$beta_y_mean)
  expect_null(p_default$beta_y_cov)
})
