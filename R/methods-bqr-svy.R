# ======================================================================
# Generic accessor functions for "bqr.svy" fits
#
# Quantiles are selected numerically through `tau`, resolved by .tau_index(),
# so users never build the internal "tau=0.500" labels themselves.
#
# Shape convention, following quantreg: a single quantile gives the ordinary
# base-R shape (a vector for coef(), a matrix for vcov()); several quantiles
# give one column per quantile, or a named list when one object cannot hold them.
# ======================================================================

#' Extract results from a fitted survey quantile regression
#'
#' Generic accessor functions for objects of class \code{"bqr.svy"} returned by
#' \code{\link{bqr.svy}}, following the conventions of the standard regression
#' methods in \pkg{stats}.
#'
#' \code{coef} returns the posterior means of the regression coefficients.
#'
#' \code{fitted} evaluates the estimated conditional quantile function
#' \eqn{\hat{Q}_\tau(y \mid x) = x^\prime \bar{\beta}(\tau)} at the observed
#' covariates.
#'
#' \code{vcov} returns the posterior covariance matrix of the coefficients,
#' computed from the posterior draws. It is a numeric matrix, stored at full
#' precision; only its \code{print} method rounds, to \code{digits} decimal
#' places.
#'
#' \code{as.matrix} returns the posterior draws themselves. They are the basis
#' for any further inference, such as posterior probabilities like
#' \eqn{P(\beta_j > 0)}, contrasts between quantiles, or passing the draws to
#' packages such as \pkg{coda} or \pkg{bayesplot}.
#'
#' \code{sigma} returns the scale parameter of the working likelihood. It is the
#' posterior mean when \code{method = "ald"} and \code{estimate_sigma = TRUE},
#' exactly \code{1} when \code{"ald"} holds it fixed, and \code{NA} for the
#' \code{"score"} and \code{"approximate"} methods, which have no scale
#' parameter.
#'
#' \code{formula} returns the formula as supplied, \code{terms} its terms
#' object, and \code{model.matrix} the design matrix the fit was built on.
#' \code{weights} returns the survey weights as supplied, without
#' normalisation. \code{update} refits the model with part of the call changed.
#'
#' These are point summaries, evaluated at the posterior mean. The posterior
#' credible intervals of the coefficients are computed by
#' \code{\link[=posterior_interval.bqr.svy]{posterior_interval}}.
#'
#' @param object,x an object of class \code{"bqr.svy"}.
#' @param tau numeric vector of fitted quantiles to extract. \code{NULL} (default)
#'   selects every quantile in the fit.
#' @param include_sigma logical; whether to keep the scale parameter column in
#'   the returned draws. Defaults to whether the model estimated it.
#' @param digits integer; number of decimal places used when printing the
#'   covariance matrix returned by \code{vcov}. Default \code{4}.
#' @param formula. a change to the model formula, in the form used by
#'   \code{\link[stats]{update}}.
#' @param evaluate logical; if \code{TRUE} (default) the updated call is
#'   evaluated and the new fit returned, otherwise the call itself is returned.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{coef} and \code{fitted} return a vector for a single quantile,
#'   otherwise a matrix with one column per quantile. \code{vcov} returns a
#'   matrix for a single quantile, otherwise a named list of matrices. Matrix
#'   arithmetic on it behaves as on any numeric matrix; \code{unclass} drops the
#'   class that governs printing.
#'   \code{as.matrix} returns the draws by parameters matrix of posterior draws
#'   for one quantile. \code{sigma} returns a scalar for a single quantile and a
#'   named vector otherwise. \code{nobs} returns the number of observations and
#'   \code{weights} the vector of survey weights.
#'
#' @seealso \code{\link{bqr.svy}}, \code{\link{bqr.svy.posterior}} for the
#'   posterior draws and credible intervals,
#'   \code{\link[=predict.bqr.svy]{predict}} for the conditional quantile at
#'   new covariate values, \code{\link{diagnostics}},
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
#' vcov(fit, tau = 0.5)
#' head(fitted(fit, tau = 0.5))
#' nobs(fit)
#' head(weights(fit))       # the survey weights, as supplied
#'
#' # the posterior draws, one row per retained iteration
#' draws <- as.matrix(fit, tau = 0.5)
#' dim(draws)
#' mean(draws[, "x1"] > 0)                       # P(beta_x1 > 0)
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
  B <- object$beta[, idx, drop = FALSE]
  X <- .model_matrix(object)[, rownames(B), drop = FALSE]
  out <- X %*% B
  if (ncol(out) == 1L) drop(out) else out
}

#' @rdname bqr.svy.methods
#' @exportS3Method stats::vcov bqr.svy
vcov.bqr.svy <- function(object, tau = NULL, ...) {
  idx <- .tau_index(object, tau)
  out <- lapply(object$quantile[idx], function(tt) {
    V <- stats::var(.draws_for_tau(object, tt, include_sigma = FALSE))
    # a plain numeric matrix, carrying a class only so that printing rounds it
    # the way the other print methods in the package do
    class(V) <- c("bqr.svy.vcov", "matrix", "array")
    V
  })
  names(out) <- .tau_labels(object)[idx]
  if (length(out) == 1L) out[[1L]] else out
}

#' @rdname bqr.svy.methods
#' @exportS3Method print bqr.svy.vcov
print.bqr.svy.vcov <- function(x, digits = 4, ...) {
  V <- unclass(x)
  fmt <- formatC(V, format = "f", digits = digits)
  dim(fmt) <- dim(V)
  dimnames(fmt) <- dimnames(V)
  print.default(fmt, quote = FALSE, right = TRUE, ...)
  invisible(x)
}

#' @rdname bqr.svy.methods
#' @exportS3Method stats::nobs bqr.svy
nobs.bqr.svy <- function(object, ...) as.integer(object$n)

#' @rdname bqr.svy.methods
#' @exportS3Method base::as.matrix bqr.svy
as.matrix.bqr.svy <- function(x, tau = NULL, include_sigma = NULL, ...) {
  .draws_for_tau(x, tau, include_sigma = include_sigma)
}

#' @rdname bqr.svy.methods
#' @exportS3Method stats::sigma bqr.svy
sigma.bqr.svy <- function(object, tau = NULL, ...) {
  idx    <- .tau_index(object, tau)
  labels <- .tau_labels(object)[idx]

  out <- if (!identical(object$method, "ald")) {
    # "score" and "approximate" have no scale parameter at all
    rep(NA_real_, length(idx))
  } else if (!isTRUE(object$estimate_sigma)) {
    # the "ald" working likelihood holds sigma at 1 when it is not estimated
    rep(1, length(idx))
  } else {
    vapply(object$draws[idx], function(m) {
      m <- as.matrix(m)
      p <- nrow(object$beta)
      if (ncol(m) > p) mean(m[, p + 1L], na.rm = TRUE) else NA_real_
    }, numeric(1))
  }

  names(out) <- labels
  if (length(out) == 1L) unname(out) else out
}


# ----------------------------------------------------------------------
# Model specification.
#
# formula(), terms(), weights() and update() already worked through the stats
# default methods, which read $formula, $terms, $weights and $call off the
# object. They are given explicit methods anyway: the default methods do not
# appear in methods(class = "bqr.svy"), so a user surveying what the class
# supports could not tell that these are available. model.matrix() is the one
# that genuinely needed a method -- model.matrix.default() tries to evaluate the
# formula in the caller and fails on the response.
# ----------------------------------------------------------------------

#' @rdname bqr.svy.methods
#' @exportS3Method stats::formula bqr.svy
formula.bqr.svy <- function(x, ...) x$formula

#' @rdname bqr.svy.methods
#' @exportS3Method stats::terms bqr.svy
terms.bqr.svy <- function(x, ...) x$terms

#' @rdname bqr.svy.methods
#' @exportS3Method stats::model.matrix bqr.svy
model.matrix.bqr.svy <- function(object, ...) .model_matrix(object)

#' @rdname bqr.svy.methods
#' @exportS3Method stats::weights bqr.svy
weights.bqr.svy <- function(object, ...) object$weights

#' @rdname bqr.svy.methods
#' @exportS3Method stats::update bqr.svy
update.bqr.svy <- function(object, formula., ..., evaluate = TRUE) NextMethod()
