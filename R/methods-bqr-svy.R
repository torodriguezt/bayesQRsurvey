# ======================================================================
# Extractor methods for "bqr.svy" fits
#
# Quantiles are selected numerically through `tau`, resolved by .tau_index(),
# so users never build the internal "tau=0.500" labels themselves.
#
# Shape convention, following quantreg: a single quantile gives the ordinary
# base-R shape (a vector for coef(), a matrix for confint()); several quantiles
# give one column per quantile, or a named list when one object cannot hold them.
# ======================================================================

#' Extract results from a fitted survey quantile regression
#'
#' Extractor methods for objects of class \code{"bqr.svy"} returned by
#' \code{\link{bqr.svy}}, following the conventions of the standard regression
#' methods in \pkg{stats}.
#'
#' \code{coef} returns the posterior means of the regression coefficients.
#'
#' \code{fitted} evaluates the estimated conditional quantile function
#' \eqn{\hat{Q}_\tau(y \mid x) = x^\prime \bar{\beta}(\tau)} at the observed
#' covariates. Note that the corresponding residuals are not centred at zero, as
#' they would be for a mean regression: because \eqn{\beta(\tau)} solves the
#' survey-weighted check-function problem, it is the \emph{weighted} share of
#' negative residuals that approaches \eqn{\tau}.
#'
#' \code{confint} returns equal-tailed credible intervals obtained from the
#' posterior draws, so its values agree with those reported by
#' \code{\link[=summary.bqr.svy]{summary}}. They describe posterior uncertainty
#' about where the conditional quantile lies, and are not prediction intervals
#' for a new observation.
#'
#' @param object an object of class \code{"bqr.svy"}.
#' @param tau numeric vector of fitted quantiles to extract. \code{NULL} (default)
#'   selects every quantile in the fit. Asking for a quantile that was not fitted
#'   is an error.
#' @param parm a specification of which parameters to report, either a vector of
#'   names or of positions. Defaults to all coefficients.
#' @param level the credible level required.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{coef} a named vector for a single quantile, otherwise a matrix
#'   with one column per quantile. \code{fitted} a vector for a single quantile,
#'   otherwise a matrix with one column per quantile. \code{confint} a matrix
#'   with the lower and upper bounds for a single quantile, otherwise a named
#'   list of such matrices.
#'
#' @seealso \code{\link{bqr.svy}}, \code{\link{diagnostics}},
#'   \code{\link[=summary.bqr.svy]{summary}}
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
#' coef(fit)                # one column per quantile
#' coef(fit, tau = 0.5)     # a single named vector
#' confint(fit, tau = 0.5)
#' head(fitted(fit, tau = 0.5))
#' }
#'
#' @name bqr.svy.methods
#' @rdname bqr.svy.methods
NULL


#' @rdname bqr.svy.methods
#' @exportS3Method stats::coef bqr.svy
coef.bqr.svy <- function(object, tau = NULL, ...) {
  idx <- .tau_index(object, tau)
  M <- object$beta[, idx, drop = FALSE]
  if (ncol(M) == 1L) stats::setNames(as.numeric(M), rownames(M)) else M
}

#' @rdname bqr.svy.methods
#' @exportS3Method stats::fitted bqr.svy
fitted.bqr.svy <- function(object, tau = NULL, ...) {
  idx <- .tau_index(object, tau)
  X <- stats::model.matrix(object$terms, object$model)
  B <- object$beta[, idx, drop = FALSE]
  out <- X[, rownames(B), drop = FALSE] %*% B
  if (ncol(out) == 1L) drop(out) else out
}

#' @rdname bqr.svy.methods
#' @exportS3Method stats::confint bqr.svy
confint.bqr.svy <- function(object, parm, level = 0.95, tau = NULL, ...) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number in (0,1).", call. = FALSE)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  idx   <- .tau_index(object, tau)

  # resolved here rather than inside the loop: missing() only works in the frame
  # where the argument is a formal
  want <- if (missing(parm)) NULL else parm

  out <- lapply(object$quantile[idx], function(tt) {
    D <- .draws_for_tau(object, tt, include_sigma = FALSE)
    if (!is.null(want)) {
      cn <- colnames(D)
      sel <- if (is.character(want)) {
        bad <- setdiff(want, cn)
        if (length(bad))
          stop("Unknown parameter(s): ", paste(bad, collapse = ", "), ".", call. = FALSE)
        want
      } else cn[want]
      D <- D[, sel, drop = FALSE]
    }
    # same computation as summary.bqr.svy, so the two agree exactly
    ci <- t(apply(D, 2, stats::quantile, probs = probs, na.rm = TRUE))
    colnames(ci) <- paste0(format(100 * probs, trim = TRUE), " %")
    ci
  })
  names(out) <- .tau_labels(object)[idx]
  if (length(out) == 1L) out[[1L]] else out
}
