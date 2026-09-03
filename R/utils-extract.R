# ======================================================================
# Shared accessors for fitted objects
#
# Every extractor method (coef, vcov, confint, fitted, predict, as.matrix,
# diagnostics, ...) and every plotting routine resolves quantiles and pulls
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
# The "ald" backend always returns a sigma column, which is a constant column of
# ones when estimate_sigma = FALSE; it is dropped unless the fit actually
# estimated it, or the caller asks for it explicitly.
.draws_for_tau <- function(x, tau = NULL, include_sigma = NULL) {
  i <- .tau_index(x, tau, single = TRUE)

  D <- as.matrix(x$draws[[i]])
  if (is.null(colnames(D)))
    colnames(D) <- paste0("V", seq_len(ncol(D)))

  keep_sigma <- if (is.null(include_sigma)) isTRUE(x$estimate_sigma) else isTRUE(include_sigma)
  if (!keep_sigma)
    D <- D[, setdiff(colnames(D), "sigma"), drop = FALSE]

  D
}
