# ======================================================================
# Generic accessor functions for "mo.bqr.svy" fits
#
# Estimation is by EM, so the object carries posterior modes and no sample from
# the posterior. The method set therefore follows what comparable EM-based
# packages provide on the fitted object itself -- point summaries only, as in
# lqmm and mclust. A fit reports point estimates and carries no quantification
# of posterior uncertainty; vcov() is defined solely so that it fails with that
# reason rather than reporting a missing method.
#
# Quantiles are resolved with the same .tau_index() and .tau_labels() helpers the
# single-output methods use; both read only $quantile, which mo.bqr.svy also has.
# ======================================================================

#' Extract results from a fitted multiple-output survey quantile regression
#'
#' Generic accessor functions for objects of class \code{"mo.bqr.svy"} returned
#' by \code{\link{mo.bqr.svy}}.
#'
#' Because the model is directional, results are indexed by quantile and by
#' direction. For a given quantile the coefficient matrix has one column per
#' direction, with rows for the covariate coefficients followed by the
#' \eqn{\gamma} coefficients on the projected responses. Quantiles are selected
#' numerically through \code{tau} and directions through \code{direction}.
#'
#' Estimation uses an EM algorithm, so a fit stores the posterior mode of the
#' coefficients rather than a sample from the posterior. Unlike a fit produced
#' by \code{\link{bqr.svy}}, it holds no posterior draws, and a covariance
#' matrix cannot be computed from point estimates alone. \code{vcov} is
#' therefore defined only to fail with a message saying so, rather than being
#' left to dispatch to the default method and return something misleading. For
#' the same reason no \code{posterior_interval} method is provided. Use
#' \code{\link{diagnostics}} to check that the algorithm converged.
#'
#' What the model estimates for a given covariate profile is a quantile region
#' of the response space, the intersection of the half-spaces fitted along each
#' direction, rather than a single value. That region is computed, drawn and
#' returned by \code{\link[=plot.mo.bqr.svy]{plot}} through its \code{xValue}
#' argument.
#'
#' \code{formula}, \code{terms} and \code{model.matrix} recover the model
#' specification, \code{weights} returns the survey weights as supplied, and
#' \code{update} refits the model with part of the call changed.
#'
#' @param object an object of class \code{"mo.bqr.svy"}.
#' @param tau numeric vector of fitted quantiles to extract. \code{NULL} (default)
#'   selects every quantile in the fit.
#' @param direction optional integer vector selecting directions. \code{NULL}
#'   (default) keeps all of them.
#' @param x an object of class \code{"mo.bqr.svy"}.
#' @param formula. a change to the model formula, in the form used by
#'   \code{\link[stats]{update}}.
#' @param evaluate logical; if \code{TRUE} (default) the updated call is
#'   evaluated and the new fit returned, otherwise the call itself is returned.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{coef} returns a matrix of coefficients, parameters by
#'   directions, for a single quantile, and a named list of such matrices
#'   otherwise. \code{sigma} returns the scale parameter for each direction,
#'   which is \code{1} throughout when \code{estimate_sigma = FALSE}. \code{nobs}
#'   returns the number of observations and \code{weights} the vector of survey
#'   weights. \code{vcov} does not return a value; it throws an error explaining
#'   that this class carries no posterior uncertainty. \code{formula},
#'   \code{terms} and \code{model.matrix} return the model specification, and
#'   \code{update} returns the refitted model.
#'
#' @seealso \code{\link{mo.bqr.svy}}, \code{\link{diagnostics}},
#'   \code{\link[=plot.mo.bqr.svy]{plot}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 250
#' d <- data.frame(x1 = runif(n, -1, 1), w = runif(n, 1, 3))
#' d$y1 <- 1 + 1.2 * d$x1 + rnorm(n)
#' d$y2 <- 0.5 - 0.7 * d$x1 + rnorm(n)
#'
#' fit <- mo.bqr.svy(cbind(y1, y2) ~ x1, weights = w, data = d,
#'                   quantile = c(0.25, 0.5), n_dir = 4)
#'
#' coef(fit, tau = 0.5)                      # parameters by directions
#' coef(fit, tau = 0.5, direction = 1:2)     # a subset of directions
#' sigma(fit, tau = 0.5)
#' nobs(fit)
#' head(weights(fit))                        # the survey weights, as supplied
#' diagnostics(fit, tau = 0.5)
#' }
#'
#' @name mo.bqr.svy.methods
#' @rdname mo.bqr.svy.methods
NULL


# Validate a direction selection against the number of directions in the fit.
.dir_index <- function(object, direction = NULL) {
  K <- object$n_dir
  if (is.null(direction)) return(seq_len(K))
  if (!is.numeric(direction) || !length(direction) || any(!is.finite(direction)))
    stop("'direction' must be a numeric vector of direction indices.", call. = FALSE)
  bad <- direction[direction < 1 | direction > K | direction != as.integer(direction)]
  if (length(bad))
    stop("Invalid direction(s): ", paste(bad, collapse = ", "),
         ". This fit has ", K, " directions.", call. = FALSE)
  as.integer(direction)
}


#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::coef mo.bqr.svy
coef.mo.bqr.svy <- function(object, tau = NULL, direction = NULL, ...) {
  idx <- .tau_index(object, tau)
  k   <- .dir_index(object, direction)

  out <- lapply(object$coefficients[idx], function(M) M[, k, drop = FALSE])
  names(out) <- .tau_labels(object)[idx]
  if (length(out) == 1L) out[[1L]] else out
}

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::sigma mo.bqr.svy
sigma.mo.bqr.svy <- function(object, tau = NULL, direction = NULL, ...) {
  idx <- .tau_index(object, tau)
  k   <- .dir_index(object, direction)

  # $sigma already holds 1 for every direction when estimate_sigma = FALSE
  out <- lapply(object$sigma[idx], function(s) s[k])
  names(out) <- .tau_labels(object)[idx]
  if (length(out) == 1L) out[[1L]] else out
}

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::nobs mo.bqr.svy
nobs.mo.bqr.svy <- function(object, ...) as.integer(object$n_obs)


# ----------------------------------------------------------------------
# Uncertainty is not available for this class.
#
# EM returns the posterior mode and no sample from the posterior, so nothing in
# the object supports a covariance matrix. vcov() exists only so that the
# generic stops with that reason instead of reporting a missing method.
#
# posterior_interval() is deliberately left unregistered: a class that cannot
# compute a credible interval should not appear to support the generic for one.
# It is the bqr.svy fits that implement it. confint() is likewise left to R:
# the package does not redefine a base generic.
# ----------------------------------------------------------------------

.mo_point_only <- function(what) {
  stop("A 'mo.bqr.svy' fit is estimated by EM and stores the posterior mode ",
       "only, so ", what, " cannot be computed from it: the object carries no ",
       "posterior sample. This class reports point estimates. Use coef() for ",
       "the estimates and diagnostics() to check that the algorithm converged.",
       call. = FALSE)
}

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::vcov mo.bqr.svy
vcov.mo.bqr.svy <- function(object, ...) .mo_point_only("a covariance matrix")


# ----------------------------------------------------------------------
# Model specification. Same five methods, same reasoning, as for "bqr.svy":
# the object stores $formula, $terms, $model, $weights and $call, so only
# model.matrix() needs real work; the rest are made explicit so that they show
# up in methods(class = "mo.bqr.svy").
#
# model.matrix() returns the covariate design matrix. The response is the
# multivariate cbind(...) on the left of the formula and is not part of it.
# weights() returns the weights as supplied; the EM works on a normalised copy.
# ----------------------------------------------------------------------

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::formula mo.bqr.svy
formula.mo.bqr.svy <- function(x, ...) x$formula

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::terms mo.bqr.svy
terms.mo.bqr.svy <- function(x, ...) x$terms

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::model.matrix mo.bqr.svy
model.matrix.mo.bqr.svy <- function(object, ...) .model_matrix(object)

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::weights mo.bqr.svy
weights.mo.bqr.svy <- function(object, ...) object$weights

#' @rdname mo.bqr.svy.methods
#' @exportS3Method stats::update mo.bqr.svy
update.mo.bqr.svy <- function(object, formula., ..., evaluate = TRUE) NextMethod()
