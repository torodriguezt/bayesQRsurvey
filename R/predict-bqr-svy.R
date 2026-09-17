# ======================================================================
# Prediction for "bqr.svy" fits.
#
# The structure follows quantreg::predict.rq, which is the reference
# implementation for this model class: no newdata returns the fitted values,
# otherwise the design matrix is rebuilt from the stored terms, xlevels and
# contrasts, and the point prediction is X %*% beta.
#
# predict.rq builds its band as X %*% t(B) over bootstrap replicates B and takes
# row quantiles. The Bayesian version is the same operation over the posterior
# draws, so the two agree in structure and differ only in where the replicates
# come from.
#
# interval = "prediction" is deliberately absent, as it is in predict.rq: a
# fitted quantile is one slice of the conditional distribution, and simulating a
# new response would require the whole of it.
# ======================================================================

#' Predicted quantiles from a fitted survey quantile regression
#'
#' Evaluates the estimated conditional quantile function
#' \eqn{\hat{Q}_\tau(y \mid x) = x^\prime \bar{\beta}(\tau)} at new covariate
#' values, optionally with a credible band.
#'
#' What is returned is the conditional quantile itself. With
#' \code{tau = 0.9}, it is the value below which 90% of the responses fall
#' among units sharing those covariates. It is a summary of the conditional
#' distribution, in the same sense that \code{\link[stats]{predict.lm}} returns
#' the conditional mean, and not a draw of a new response.
#'
#' \code{interval = "credible"} adds a band obtained from the posterior draws.
#' Each draw of \eqn{\beta(\tau)} gives a draw of \eqn{x^\prime\beta(\tau)}, and
#' \code{lower} and \code{higher} are the percentiles of those values at the
#' requested \code{level}. It describes uncertainty
#' about where the conditional quantile lies, and its width reflects only that;
#' individual responses scatter much more widely around it.
#'
#' There is no \code{interval = "prediction"}, following
#' \code{\link[quantreg]{predict.rq}}. A prediction interval for a new
#' observation requires the whole conditional distribution of the response,
#' whereas a quantile regression estimates one slice of it. For
#' \code{method = "ald"} the asymmetric Laplace working likelihood is formally a
#' complete distribution, but its shape away from \eqn{\tau} is an artefact of
#' the estimation device rather than a claim about the data.
#'
#' When several quantiles are requested, each was fitted separately and nothing
#' constrains them to be ordered, so predicted quantiles may cross at some
#' covariate values.
#'
#' @param object an object of class \code{"bqr.svy"}.
#' @param newdata an optional data frame of covariate values. If omitted, the
#'   fitted values are returned.
#' @param tau numeric vector of fitted quantiles to predict. \code{NULL}
#'   (default) uses every quantile in the fit.
#' @param interval \code{"none"} (default) for point predictions, or
#'   \code{"credible"} to add a band from the posterior draws.
#' @param level the credible level of the band.
#' @param na.action how to handle missing values in \code{newdata}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return With \code{interval = "none"}, a vector for a single quantile and a
#'   matrix with one column per quantile otherwise. With
#'   \code{interval = "credible"}, a matrix with columns \code{fit},
#'   \code{lower} and \code{higher} for a single quantile, and a named list of
#'   such matrices otherwise.
#'
#' @seealso \code{\link{bqr.svy.methods}},
#'   \code{\link[=posterior_interval.bqr.svy]{posterior_interval}},
#'   \code{\link[quantreg]{predict.rq}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 300
#' d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), w = runif(n, 1, 3))
#' d$y <- 2 + 1.5 * d$x1 - 0.8 * d$x2 + rnorm(n)
#'
#' fit <- bqr.svy(y ~ x1 + x2, weights = w, data = d,
#'                quantile = c(0.1, 0.5, 0.9), niter = 2000, burnin = 500)
#'
#' nd <- data.frame(x1 = c(-0.5, 0, 0.5), x2 = 0)
#'
#' predict(fit, nd, tau = 0.9)
#' predict(fit, nd, tau = 0.9, interval = "credible")
#' predict(fit, nd)                      # one column per fitted quantile
#' }
#'
#' @exportS3Method stats::predict bqr.svy
predict.bqr.svy <- function(object, newdata, tau = NULL,
                            interval = c("none", "credible"),
                            level = 0.95, na.action = stats::na.pass, ...) {
  interval <- match.arg(interval)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number in (0,1).", call. = FALSE)

  idx    <- .tau_index(object, tau)
  labels <- .tau_labels(object)[idx]
  B      <- object$beta

  if (missing(newdata) || is.null(newdata)) {
    X <- .model_matrix(object)
  } else {
    Terms <- stats::delete.response(object$terms)
    mf <- stats::model.frame(Terms, newdata, na.action = na.action,
                             xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      stats::.checkMFClasses(cl, mf)
    X <- stats::model.matrix(Terms, mf, contrasts.arg = object$contrasts)
  }

  missing_cols <- setdiff(rownames(B), colnames(X))
  if (length(missing_cols))
    stop("The design matrix of 'newdata' is missing: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  X <- X[, rownames(B), drop = FALSE]

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)

  one_tau <- function(i) {
    fit <- drop(X %*% B[, i])
    if (interval == "none") return(fit)

    # the same operation predict.rq performs over bootstrap replicates
    D  <- .draws_for_tau(object, object$quantile[i], include_sigma = FALSE)
    XB <- X %*% t(D[, rownames(B), drop = FALSE])
    band <- t(apply(XB, 1, stats::quantile, probs = probs, na.rm = TRUE))
    out <- cbind(fit, band)
    colnames(out) <- c("fit", "lower", "higher")
    out
  }

  res <- lapply(idx, one_tau)
  names(res) <- labels
  if (length(res) == 1L) return(res[[1L]])
  if (interval == "none") {
    out <- do.call(cbind, res)
    colnames(out) <- labels
    return(out)
  }
  res
}
