# ======================================================================
# Credible intervals, through the rstantools generic that brms and rstanarm
# build on.
#
# Naming the method posterior_interval() rather than confint() keeps the
# Bayesian reading explicit: what is returned is a credible interval, not a
# confidence interval.
# ======================================================================

#' Credible intervals from a fitted survey quantile regression
#'
#' Computes equal-tailed credible intervals for the regression coefficients of
#' a fitted \code{"bqr.svy"} model from its posterior draws.
#'
#' The values agree with those reported by
#' \code{\link[=summary.bqr.svy]{summary}}. They describe posterior uncertainty
#' about where the conditional quantile lies, and are not prediction intervals
#' for a new observation.
#'
#' \code{posterior_interval} is the generic defined in \pkg{rstantools} and
#' used by \pkg{brms} and \pkg{rstanarm}, so the same call works across these
#' packages.
#'
#' @param object an object of class \code{"bqr.svy"}.
#' @param tau a single fitted quantile. Defaults to the only one in the fit, and
#'   is required when several were fitted.
#' @param prob the credible level. Default \code{0.95}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A matrix with one row per coefficient and columns holding the lower
#'   and upper bounds.
#'
#' @seealso \code{\link{bqr.svy.methods}}, \code{\link{diagnostics}},
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
#'                quantile = c(0.25, 0.5), niter = 2000, burnin = 500)
#'
#' posterior_interval(fit, tau = 0.5)
#' posterior_interval(fit, tau = 0.5, prob = 0.9)
#' }
#'
#' @name bqr.svy.posterior
#' @rdname bqr.svy.posterior
NULL


#' @rdname bqr.svy.posterior
#' @importFrom rstantools posterior_interval
#' @exportS3Method rstantools::posterior_interval bqr.svy
posterior_interval.bqr.svy <- function(object, prob = 0.95, tau = NULL, ...) {
  if (!is.numeric(prob) || length(prob) != 1L || prob <= 0 || prob >= 1)
    stop("'prob' must be a single number in (0,1).", call. = FALSE)

  probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  D <- .draws_for_tau(object, tau, include_sigma = FALSE)

  # the same stats::quantile path summary.bqr.svy uses, so the two agree exactly
  ci <- t(apply(D, 2, stats::quantile, probs = probs, na.rm = TRUE))
  colnames(ci) <- paste0(format(100 * probs, trim = TRUE), "%")
  ci
}


# Re-export the generics so they are callable without attaching rstantools.

#' @importFrom rstantools posterior_interval
#' @export
rstantools::posterior_interval
