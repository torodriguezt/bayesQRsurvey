# ======================================================================
# Convergence diagnostics
#
# Public replacement for reaching into fit$diagnosis[["tau=0.500"]].
# ======================================================================

#' Convergence diagnostics for fitted survey quantile regressions
#'
#' Extracts the convergence diagnostics stored by \code{\link{bqr.svy}} and
#' \code{\link{mo.bqr.svy}}, selecting quantiles numerically rather than by
#' internal label.
#'
#' For \code{"bqr.svy"} fits the diagnostics are those of the MCMC sampler: the
#' potential scale reduction factor \eqn{\hat{R}} together with the bulk and tail
#' effective sample sizes, computed with the \pkg{posterior} package. As a rule of
#' thumb \eqn{\hat{R} < 1.01} and effective sample sizes above 400 indicate that
#' the chain has converged and carries enough information for stable posterior
#' summaries.
#'
#' For \code{"mo.bqr.svy"} fits estimation is by an EM algorithm rather than MCMC,
#' so the diagnostics report, for each quantile and direction, the number of
#' iterations taken and whether the tolerance was met before \code{max_iter}.
#'
#' @param object a fitted object of class \code{"bqr.svy"} or \code{"mo.bqr.svy"}.
#' @param tau numeric vector of fitted quantiles. \code{NULL} (default) selects
#'   every quantile in the fit.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A data frame of diagnostics for a single quantile, or a named list of
#'   such data frames when several quantiles are selected.
#'
#' @seealso \code{\link{bqr.svy}}, \code{\link{mo.bqr.svy}},
#'   \code{\link{bqr.svy.methods}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 300
#' d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), w = runif(n, 1, 3))
#' d$y <- 2 + 1.5 * d$x1 - 0.8 * d$x2 + rnorm(n)
#'
#' fit <- bqr.svy(y ~ x1 + x2, weights = w, data = d,
#'                quantile = c(0.25, 0.5, 0.75), niter = 2000, burnin = 500)
#'
#' diagnostics(fit, tau = 0.5)
#' diagnostics(fit)              # all fitted quantiles
#' }
#'
#' @export
diagnostics <- function(object, ...) UseMethod("diagnostics")

#' @rdname diagnostics
#' @exportS3Method diagnostics bqr.svy
diagnostics.bqr.svy <- function(object, tau = NULL, ...) {
  idx <- .tau_index(object, tau)
  out <- object$diagnosis[idx]
  names(out) <- .tau_labels(object)[idx]
  if (length(out) == 1L) out[[1L]] else out
}

#' @rdname diagnostics
#' @exportS3Method diagnostics mo.bqr.svy
diagnostics.mo.bqr.svy <- function(object, tau = NULL, ...) {
  taus <- object$quantile
  idx <- if (is.null(tau)) {
    seq_along(taus)
  } else {
    vapply(tau, function(tt) {
      d <- abs(taus - tt)
      j <- which.min(d)
      if (d[j] > 1e-8)
        stop("Quantile ", format(tt), " was not fitted. Available quantiles: ",
             paste(format(taus), collapse = ", "), ".", call. = FALSE)
      j
    }, integer(1))
  }

  labels <- paste0("tau=", formatC(taus, format = "f", digits = 3))
  out <- lapply(idx, function(i) {
    dirs <- object$fit[[i]]$directions
    res <- data.frame(
      direction = names(dirs),
      iter      = vapply(dirs, function(z) as.numeric(z$iter), numeric(1)),
      converged = vapply(dirs, function(z) isTRUE(as.logical(z$converged)), logical(1)),
      sigma     = vapply(dirs, function(z) as.numeric(z$sigma)[1], numeric(1)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    class(res) <- c("bqr.svy.diagnostics", "data.frame")
    res
  })
  names(out) <- labels[idx]
  if (length(out) == 1L) out[[1L]] else out
}

#' @exportS3Method print bqr.svy.diagnostics
print.bqr.svy.diagnostics <- function(x, digits = 3, ...) {
  df <- as.data.frame(x)
  if ("rhat" %in% names(df))
    df$rhat <- formatC(df$rhat, format = "f", digits = digits)
  for (nm in intersect(c("ess_bulk", "ess_tail"), names(df)))
    df[[nm]] <- round(df[[nm]])
  if ("sigma" %in% names(df))
    df$sigma <- formatC(df$sigma, format = "f", digits = digits)
  print.data.frame(df, ...)
  invisible(x)
}
