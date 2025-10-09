# =============================================================================
# Prior Interface for bayesQRsurvey
# =============================================================================

if (!exists("%||%"))
  `%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# Small helper: expand scalar/vector to diagonal p x p
.as_diag <- function(x, p, what) {
  if (is.null(x)) return(NULL)
  if (is.numeric(x) && length(x) == 1L) return(diag(x, p))
  if (is.numeric(x) && is.null(dim(x)) && length(x) == p) return(diag(x, p))
  if (is.matrix(x)) {
    if (all(dim(x) == c(p, p))) return(x)
    stop(what, " must be a ", p, "x", p, " matrix.", call. = FALSE)
  }
  stop(what, " must be scalar, length-", p, " vector, or ", p, "x", p, " matrix.", call. = FALSE)
}

#' Create prior for Bayesian quantile regression models for complex survey data
#'
#' \code{prior} creates prior distributions for both single (\code{bqr.svy}) and multiple-output
#' (\code{mo.bqr.svy}) Bayesian quantile regression models for complex survey data.
#'
#' @param beta_x_mean vector of prior means for the regression coefficients. (default = NULL).
#' @param beta_x_cov prior covariance matrix for the regression coefficients. (default = NULL).
#' @param sigma_shape shape parameter for inverse Gamma prior for \eqn{\sigma^2}. (default = 0.001).
#' @param sigma_rate rate parameter for inverse Gamma prior for \eqn{\sigma^2}. (default = 0.001).
#' @param beta_y_mean prior means for the coefficients related to the variables that emerge from the product between the orthogonal basis and the outputs
#' (default = NULL).
#' @param beta_y_cov prior covariance matrix for the coefficients related to the variables that emerge from the product between the orthogonal basis and the outputs.
#'  (default = NULL).
#'
#' @details
#' The function \code{prior} builds prior distributions for the three methods implemented in the function
#' \code{bqr.svy} and for the multiple-output quantile regression implemented in the function \code{mo.bqr.svy}.
#' Every nonspecified prior parameter will get the default value.
#'
#' \itemize{
#'   \item \code{method = "ald"} in function \code{bqr.svy} allow the specification of hyperparameters
#'         \code{beta_x_mean}, \code{beta_x_cov}, \code{sigma_shape}, and \code{sigma_rate}.
#'   \item \code{method = "score"} in function \code{bqr.svy} allow the specification of hyperparameters
#'         \code{beta_x_mean} and \code{beta_x_cov}.
#'   \item \code{method = "approximate"} in function \code{bqr.svy} allow the specification of hyperparameters
#'         \code{beta_x_mean} and \code{beta_x_cov}.
#'   \item In function \code{mo.bqr.svy}, the specification of hyperparameters \code{beta_x_mean},\code{beta_x_cov},
#'         \code{sigma_shape}, \code{sigma_rate}, \code{beta_y_mean}, and \code{beta_y_cov} are allowed.
#' }
#'
#' @return An object of class \code{"prior"}.
#' @examples
#'
#' #Simulate data
#' set.seed(123)
#' n  <- 200
#' x1 <- rnorm(n, 0, 1)
#' x2 <- runif(n, -1, 1)
#' w  <- runif(n, 0.5, 2)   # survey weights
#'
#' y1 <- 2 + 1.5*x1 - 0.8*x2 + rnorm(n, 0, 1)
#' y2 <- 1 + 0.5*x1 - 0.2*x2 + rnorm(n, 0, 1)
#'
#' data <- data.frame(y1 = y1, y2 = y2, x1 = x1, x2 = x2, w = w)
#'
#'
#' # Define a general informative prior
#' prior_general <- prior(
#'   beta_x_mean = c(2, 1.5, -0.8),
#'   beta_x_cov  = diag(c(0.25, 0.25, 0.25)),
#'   sigma_shape = 3,
#'   sigma_rate  = 2,
#'   beta_y_mean = 1,
#'   beta_y_cov  = 0.25
#' )
#'
#' #Estimate the model parameters with informative prior
#'
#'
#' fit_ald <- bqr.svy(y1 ~ x1 + x2, weights = w, data = data,
#'                    prior = prior_general, method = "ald")
#'
#' fit_scr <- bqr.svy(y1 ~ x1 + x2, weights = w, data = data,
#'                    prior = prior_general, method = "score")
#'
#' fit_apx <- bqr.svy(y1 ~ x1 + x2, weights = w, data = data,
#'                    prior = prior_general, method = "approximate")
#'
#' # Multiple-output method
#' fit_mo <- mo.bqr.svy(cbind(y1, y2) ~ x1 + x2, weights = w,
#'                      data = data, prior = prior_general, n_dir = 10)
#'
#' plot(fit_ald, type = "trace", which = "x1", tau = 0.5)
#' 
#' plot(fit_ald, type = "trace", which = "x2", tau = 0.5)
#' 
#' print(fit_mo)
#'
#' @seealso \code{\link{bqr.svy}}, \code{\link{mo.bqr.svy}},
#'   \code{\link{summary}}
#' @export
prior <- function(
    beta_x_mean = NULL,
    beta_x_cov  = NULL,
    sigma_shape = 0.001,
    sigma_rate  = 0.001,
    beta_y_mean = NULL,
    beta_y_cov  = NULL
) {
  if (!is.null(sigma_shape) && (!is.numeric(sigma_shape) || length(sigma_shape)!=1L || sigma_shape<=0))
    stop("'sigma_shape' must be a positive scalar.", call. = FALSE)
  if (!is.null(sigma_rate) && (!is.numeric(sigma_rate) || length(sigma_rate)!=1L || sigma_rate<=0))
    stop("'sigma_rate' must be a positive scalar.", call. = FALSE)

  mc <- match.call()
  user_defined_sigma_prior <- ("sigma_shape" %in% names(mc)) || ("sigma_rate" %in% names(mc))

  structure(
    list(
      beta_x_mean = beta_x_mean,
      beta_x_cov  = beta_x_cov,
      sigma_shape = sigma_shape,
      sigma_rate  = sigma_rate,
      beta_y_mean = beta_y_mean,
      beta_y_cov  = beta_y_cov
    ),
    class = "prior",
    user_defined_sigma_prior = user_defined_sigma_prior
  )
}


#' @export
print.prior <- function(x, digits = 3, ...) {
  cat("Prior")

  # beta_x_mean
  if (is.null(x$beta_x_mean)) {
    cat("  beta_x_mean: diffuse (all zeros)\n")
  } else {
    cat("  beta_x_mean: ", paste(round(x$beta_x_mean, digits), collapse=" "), "\n", sep = "")
  }

  # beta_x_cov
  if (is.null(x$beta_x_cov)) {
    cat("  beta_x_cov:  diffuse (large diagonal: 1e6)\n")
  } else {
    cat("  beta_x_cov:\n")
    print(round(x$beta_x_cov, digits))
  }

  # sigma prior
  if (is.null(x$sigma_shape) || is.null(x$sigma_rate)) {
    cat("  sigma IG:    diffuse (shape=0.001, rate=0.001)\n")
  } else {
    cat("  sigma IG:    shape=", round(x$sigma_shape, digits),
        ", rate=", round(x$sigma_rate, digits), "\n", sep = "")
  }

  # beta_y_mean
  if (is.null(x$beta_y_mean)) {
    cat("  beta_y_mean: diffuse (all zeros)\n")
  } else {
    cat("  beta_y_mean: ", paste(round(x$beta_y_mean, digits), collapse=" "), "\n", sep = "")
  }

  # beta_y_cov
  if (is.null(x$beta_y_cov)) {
    cat("  beta_y_cov:  diffuse (large diagonal: 1e6)\n")
  } else {
    cat("  beta_y_cov:\n")
    print(round(x$beta_y_cov, digits))
  }

  invisible(x)
}



#' @keywords internal
#' @noRd
as_bqr_prior <- function(x, ...) UseMethod("as_bqr_prior")

#' @keywords internal
#' @noRd
as_mo_bqr_prior <- function(x, ...) UseMethod("as_mo_bqr_prior")

#' @keywords internal
#' @noRd
#' @exportS3Method as_bqr_prior prior
as_bqr_prior.prior <- function(x, p, names_x = NULL, ...) {
  bmean <- if (is.null(x$beta_x_mean)) rep(0, p) else if (length(x$beta_x_mean)==1L) rep(x$beta_x_mean, p) else x$beta_x_mean
  if (length(bmean)!=p) stop("beta_x_mean must have length p (", p, ").")
  bcov <- if (is.null(x$beta_x_cov)) diag(1e6, p) else .as_diag(x$beta_x_cov, p, "beta_x_cov")
  if (!is.null(names_x)) { names(bmean) <- names_x; dimnames(bcov) <- list(names_x,names_x) }
  c0 <- x$sigma_shape %||% 0.001
  C0 <- x$sigma_rate  %||% 0.001
  structure(list(b0=bmean, B0=bcov, c0=c0, C0=C0), class="bqr_prior")
}

#' @keywords internal
#' @noRd
#' @exportS3Method as_mo_bqr_prior prior
as_mo_bqr_prior.prior <- function(x, p, d, names_x=NULL, names_y=NULL, ...) {
  bmean <- if (is.null(x$beta_x_mean)) rep(0, p) else if (length(x$beta_x_mean)==1L) rep(x$beta_x_mean, p) else x$beta_x_mean
  if (length(bmean)!=p) stop("beta_x_mean must have length p (", p, ").")
  bcov <- if (is.null(x$beta_x_cov)) diag(1e6, p) else .as_diag(x$beta_x_cov, p, "beta_x_cov")

  q <- max(d-1L,0L)
  ymean <- if (q==0) NULL else if (is.null(x$beta_y_mean)) rep(0,q) else if (length(x$beta_y_mean)==1L) rep(x$beta_y_mean,q) else x$beta_y_mean
  if (q>0 && length(ymean)!=q) stop("beta_y_mean must have length d-1 (", q, ").")
  ycov <- if (q==0) NULL else if (is.null(x$beta_y_cov)) diag(1e6,q) else .as_diag(x$beta_y_cov,q,"beta_y_cov")

  if (!is.null(names_x)) { names(bmean) <- names_x; dimnames(bcov) <- list(names_x,names_x) }
  if (!is.null(ymean) && !is.null(names_y)) { names(ymean) <- names_y; if (!is.null(ycov)) dimnames(ycov) <- list(names_y,names_y) }

  structure(list(beta_mean=bmean, beta_cov=bcov,
                 sigma_shape=x$sigma_shape %||% 0.001,
                 sigma_rate =x$sigma_rate  %||% 0.001,
                 beta_star_mean=ymean, beta_star_cov=ycov),
            class="mo_bqr_prior")
}
