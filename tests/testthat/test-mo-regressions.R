# An asymmetric, tilted region with nonzero orthogonal coefficients:
# -3 <= y1 - 0.2*y2 <= 0 and 8 <= y2 - 0.3*y1 <= 12.
# Its exact coefficients isolate plotting from EM convergence.
.mo_region_fixture <- function() {
  dat <- data.frame(y1 = c(-2, 4, 0, 1), y2 = c(6, 14, 9, 11))
  U <- rbind(c(1, 0, -1, 0), c(0, 1, 0, -1))
  G <- lapply(seq_len(ncol(U)), function(k)
    matrix(c(-U[2, k], U[1, k]), ncol = 1L))
  B <- rbind(
    "(Intercept)" = c(-3, 8, 0, -12),
    gamma_1 = c(0.2, -0.3, 0.2, -0.3)
  )
  colnames(B) <- paste0("dir_", 1:4)
  structure(list(
    model = model.frame(cbind(y1, y2) ~ 1, data = dat),
    n_vars = 1L, n_dir = 4L, response_dim = 2L,
    U = U, Gamma_list = G, quantile = 0.1,
    coefficients = setNames(list(B), "tau=0.100"),
    fit = setNames(list(list()), "tau=0.100")
  ), class = "mo.bqr.svy")
}

test_that("reversing responses transposes the full fitted geometry", {
  obj <- .mo_region_fixture()
  old_y1 <- seq(-2, 4, length.out = 41)
  old_y2 <- seq(6, 14, length.out = 41)

  # Independent analytic membership for the known parallelogram.
  grid <- expand.grid(y1 = old_y1, y2 = old_y2)
  z1 <- grid$y1 - 0.2 * grid$y2
  z2 <- grid$y2 - 0.3 * grid$y1
  inside <- grid[z1 >= -3 & z1 <= 0 & z2 >= 8 & z2 <= 12, ]

  check_slices <- function(result, first, second) {
    axis <- sort(unique(first))
    lo <- vapply(axis, function(v) min(second[first == v]), numeric(1))
    hi <- vapply(axis, function(v) max(second[first == v]), numeric(1))
    expect_equal(result$data$y1, axis)
    expect_equal(result$data$min, unname(lo))
    expect_equal(result$data$max, unname(hi))
  }

  original <- plot(obj, ngridpoints = 41,
                    range_y = rbind(c(-2, 4), c(6, 14)))
  reversed <- plot(obj, response = c("y2", "y1"), ngridpoints = 41,
                    range_y = rbind(c(6, 14), c(-2, 4)))
  expect_s3_class(original$plot, "ggplot")
  expect_s3_class(reversed$plot, "ggplot")
  check_slices(original, inside$y1, inside$y2)
  check_slices(reversed, inside$y2, inside$y1)
  expect_equal(reversed$plot$labels$x, "y2")
  expect_equal(reversed$plot$labels$y, "y1")
})

test_that("response selection cannot silently change the fitted responses", {
  obj <- .mo_region_fixture()
  for (response in list(c("y1", "y1"), c("other", "y2"), c(NA_character_, "y2")))
    expect_error(plot(obj, response = response), "two fitted response names")
})

test_that("higher-dimensional fitting remains supported but plotting is bivariate", {
  set.seed(908)
  dat <- data.frame(y1 = rnorm(60), y2 = 10 + rnorm(60), y3 = rnorm(60))
  fit <- mo.bqr.svy(cbind(y1, y2, y3) ~ 1, data = dat,
                    quantile = 0.1, U = diag(3), max_iter = 500)
  expect_s3_class(fit, "mo.bqr.svy")
  expect_equal(fit$response_dim, 3L)
  expect_true(all(is.finite(coef(fit))))
  expect_error(plot(fit), "exactly two response variables")
  expect_error(plot(fit, response = c("y1", "y2")), "exactly two response variables")
})

test_that("missing model values are rejected before weights are aligned", {
  dat <- data.frame(y1 = 1:8, y2 = 11:18, x = seq_len(8), w = seq_len(8))
  dat$x[2] <- NA_real_
  expect_error(mo.bqr.svy(cbind(y1, y2) ~ x, data = dat, weights = w),
               "Data contains missing values")
  expect_error(mo.bqr.svy(cbind(y1, y2) ~ x, data = dat),
               "Data contains missing values")
  dat$x[2] <- 2
  dat$y2[3] <- NA_real_
  expect_error(mo.bqr.svy(cbind(y1, y2) ~ x, data = dat, weights = w),
               "Data contains missing values")
})

test_that("multiple-output fits reject unsupported offsets explicitly", {
  dat <- data.frame(y1 = 1:8, y2 = 11:18, x = seq_len(8), z = rep(10, 8))
  expect_error(mo.bqr.svy(cbind(y1, y2) ~ x + offset(z), data = dat),
               "offset terms are not supported")
})
