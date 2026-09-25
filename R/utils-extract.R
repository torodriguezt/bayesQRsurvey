# ======================================================================
# Shared accessors for fitted objects
#
# Every accessor method (coef, fitted, vcov, posterior_interval, diagnostics)
# and every plotting routine resolves quantiles and pulls
# posterior draws through the two helpers below, so that the "which tau is this"
# and "does this fit carry a sigma column" conventions live in exactly one place.
# ======================================================================

# Canonical per-quantile labels used to name the $draws and $diagnosis lists.
.tau_labels <- function(x) {
  paste0("tau=", formatC(x$quantile, format = "f", digits = 3))
}

# Resolve user-supplied quantiles to positions in x$quantile.
#
# `tau = NULL` means "all fitted quantiles". Matching is numeric rather than by
# label, so callers never have to reconstruct the "tau=0.500" strings, and an
# unfitted quantile is an error instead of silently snapping to the nearest one.
.tau_index <- function(x, tau = NULL, single = FALSE) {
  taus <- x$quantile

  if (is.null(tau)) {
    if (single && length(taus) != 1L)
      stop("This fit contains ", length(taus), " quantiles; specify 'tau', one of: ",
           paste(format(taus), collapse = ", "), ".", call. = FALSE)
    return(seq_along(taus))
  }

  if (!is.numeric(tau) || !length(tau) || any(!is.finite(tau)))
    stop("'tau' must be a numeric vector of fitted quantiles.", call. = FALSE)
  if (single && length(tau) != 1L)
    stop("'tau' must be a single quantile here, but ", length(tau), " were given.",
         call. = FALSE)

  vapply(tau, function(tt) {
    d <- abs(taus - tt)
    j <- which.min(d)
    if (d[j] > 1e-8)
      stop("Quantile ", format(tt), " was not fitted. Available quantiles: ",
           paste(format(taus), collapse = ", "), ".", call. = FALSE)
    j
  }, integer(1))
}

# Posterior draws for one quantile, as a matrix with named columns.
#
# Coefficients occupy the first p columns. The "ald" backend appends its scale
# in column p + 1, which is constant when estimate_sigma = FALSE. A coefficient
# may itself be named sigma, so the scale must never be selected by that name.
.sigma_draw_name <- function(coef_names) {
  make.unique(c(coef_names, "sigma"))[length(coef_names) + 1L]
}

.draws_for_tau <- function(x, tau = NULL, include_sigma = NULL) {
  i <- .tau_index(x, tau, single = TRUE)

  D <- as.matrix(x$draws[[i]])
  if (is.null(colnames(D)))
    colnames(D) <- paste0("V", seq_len(ncol(D)))
  p <- nrow(x$beta)
  if (ncol(D) < p)
    stop("The draws do not contain all fitted coefficients.", call. = FALSE)
  colnames(D)[seq_len(p)] <- rownames(x$beta)

  has_sigma <- identical(x$method, "ald") && ncol(D) > p
  if (has_sigma)
    colnames(D)[p + 1L] <- .sigma_draw_name(rownames(x$beta))

  keep_sigma <- if (is.null(include_sigma)) isTRUE(x$estimate_sigma) else isTRUE(include_sigma)
  keep <- seq_len(p)
  if (has_sigma && keep_sigma) keep <- c(keep, p + 1L)
  D[, keep, drop = FALSE]
}

# Design matrix of a fit, rebuilt from the stored terms and model frame.
#
# Both classes keep $terms and $model, so one helper serves them. Rebuilding
# rather than storing X keeps the fitted object smaller and guarantees the
# matrix always matches the terms actually used.
.model_matrix <- function(x) {
  stats::model.matrix(x$terms, x$model, contrasts.arg = x$contrasts)
}
