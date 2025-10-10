if (!exists("%||%"))
  `%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

.user_defined_sigma <- function(pr) {
  if (is.null(pr)) return(FALSE)

  uds <- attr(pr, "user_defined_sigma_prior")
  if (!is.null(uds)) return(isTRUE(uds))

  DEFAULT_SHAPE <- 0.001
  DEFAULT_RATE  <- 0.001

  if (inherits(pr, "prior")) {
    has_shape <- !is.null(pr$sigma_shape)
    has_rate  <- !is.null(pr$sigma_rate)

    if (!has_shape && !has_rate) return(FALSE)
    if (has_shape && has_rate &&
        isTRUE(all(is.finite(c(pr$sigma_shape, pr$sigma_rate))))) {
      if (identical(pr$sigma_shape, DEFAULT_SHAPE) &&
          identical(pr$sigma_rate,  DEFAULT_RATE)) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    }
    return(TRUE)
  }

  if (inherits(pr, "bqr_prior")) {
    has_c0 <- !is.null(pr$c0)
    has_C0 <- !is.null(pr$C0)
    if (!has_c0 && !has_C0) return(FALSE)
    if (has_c0 && has_C0 &&
        isTRUE(all(is.finite(c(pr$c0, pr$C0))))) {
      if (identical(pr$c0, DEFAULT_SHAPE) && identical(pr$C0, DEFAULT_RATE)) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    }
    return(TRUE)
  }

  FALSE
}

# ==== MODEL FITTER ============================================================

#' Bayesian quantile regression for complex survey data
#'
#' \code{bqr.svy} implements Bayesian methods for estimating quantile regression models
#' for complex survey data analysis regarding single (univariate) outputs. To
#' improve computational efficiency, the Markov Chain Monte Carlo (MCMC) algorithms
#' are implemented in C++.
#'
#' @param formula a symbolic description of the model to be fit.
#' @param weights an optional numerical vector containing the survey weights. If \code{NULL}, equal weights are used.
#' @param data an optional data frame containing the variables in the model.
#' @param quantile numerical scalar or vector containing quantile(s) of interest (default=0.5).
#' @param method one of \code{"ald"}, \code{"score"} and \code{"approximate"} (default=\code{"ald"}).
#' @param prior a \code{bqr_prior} object of class "prior". If omitted, a vague prior is assumed (see \code{\link{prior}}).
#' @param niter number of MCMC draws.
#' @param burnin number of initial MCMC draws to be discarded.(default = 0)
#' @param thin thinning parameter, i.e., keep every keepth draw (default=1).
#' @param verbose logical flag indicating whether to print progress messages (default=TRUE).
#' @param estimate_sigma logical flag indicating whether to estimate the scale parameter
#' when method = "ald" (default=FALSE and \eqn{\sigma^2} is set to 1)
#'
#' @details
#' The bqr.svy function can estimate three types of models, where the quantile regression
#' coefficients are defined at the super-population level, and their estimators are built upon
#' the survey weights.
#' \itemize{
#'   \item \code{"ald"} – The asymmetric Laplace distribution as working likelihood.
#'   \item \code{"score"} – A score based likelihood function.
#'   \item \code{"approximate"} – A pseudolikelihood function based on a Gaussian approximation.
#' }
#'
#' @return An object of class \code{"bqr.svy"}, containing:
#' \item{beta}{Posterior mean estimates of regression coefficients.}
#' \item{draws}{Posterior draws from the MCMC sampler.}
#' \item{accept_rate}{Average acceptance rate (if available).}
#' \item{warmup, thin}{MCMC control parameters used during sampling.}
#' \item{quantile}{The quantile(s) fitted.}
#' \item{prior}{Prior specification used.}
#' \item{formula, terms, model}{Model specification details.}
#' \item{runtime}{Elapsed runtime in seconds.}
#' \item{method}{Estimation method}
#' \item{estimate_sigma}{Logical flag indicating whether the scale parameter
#'   \eqn{\sigma^2} was estimated (\code{TRUE}) or fixed at 1 (\code{FALSE}).}
#'
#' @references
#' Nascimento, M. L. & Gon\enc{c}{c}alves, K. C. M. (2024). Bayesian Quantile Regression Models
#'   for Complex Survey Data Under Informative Sampling. *Journal of Survey Statistics and Methodology*,
#'   12(4), 1105–1130.
#'
#' @examples
#' # Generate population data
#' set.seed(123)
#' N    <- 10000
#' x1_p <- runif(N, -1, 1)
#' x2_p <- runif(N, -1, 1)
#' y_p  <- 2 + 1.5 * x1_p - 0.8 * x2_p + rnorm(N)
#'
#' # Generate sample data
#' n <- 500
#' z_aux <- rnorm(N, mean = 1 + y_p, sd=.5)
#' p_aux <- 1 / (1 + exp(2.5 - 0.5 * z_aux))
#' s_ind <- sample(1:N, n, replace = FALSE, prob = p_aux)
#' y_s   <- y_p[s_ind]
#' x1_s  <- x1_p[s_ind]
#' x2_s  <- x2_p[s_ind]
#' w     <- 1 / p_aux[s_ind]
#' data  <- data.frame(y = y_s, x1 = x1_s, x2 = x2_s, w = w)
#'
#' # Basic usage with default method ('ald') and priors (vague)
#' fit1 <- bqr.svy(y ~ x1 + x2, weights = w, data = data)
#'
#' # Specify informative priors
#' prior<- prior(
#'  beta_x_mean = c(2, 1.5, -0.8),
#'  beta_x_cov = diag(c(0.25, 0.25, 0.25)),
#'  sigma_shape = 1,
#'  sigma_rate = 1
#')
#' fit2 <- bqr.svy(y ~ x1 + x2, weights = w, data = data, prior = prior)
#'
#' # Specify different methods
#' fit_score  <- bqr.svy(y ~ x1 + x2, weights = w, data = data, method = "score")
#' fit_approx <- bqr.svy(y ~ x1 + x2, weights = w, data = data, method = "approximate")
#'
#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @export
bqr.svy <- function(formula,
                    weights  = NULL,
                    data     = NULL,
                    quantile = 0.5,
                    method   = c("ald", "score", "approximate"),
                    prior    = NULL,
                    niter    = 50000,
                    burnin   = 0,
                    thin     = 1,
                    verbose  = TRUE,
                    estimate_sigma = FALSE) {

  tic    <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  cl <- match.call()

  ## --- Aviso sobre estimate_sigma SOLO si el usuario lo pasó y el método no es 'ald'
  if (method != "ald" && !missing(estimate_sigma)) {
    warning("'estimate_sigma' only applies to the 'ald' method and will be ignored", call. = FALSE)
  }

  # --- Quantiles ---
  if (!is.numeric(quantile) || any(!is.finite(quantile)))
    stop("'quantile' must be numeric and finite.", call. = FALSE)
  if (any(quantile <= 0 | quantile >= 1))
    stop("All elements of 'quantile' must be in (0,1).", call. = FALSE)
  taus <- sort(unique(as.numeric(quantile)))
  if (length(taus) < length(quantile))
    warning("Duplicated quantiles were provided; using unique sorted values.")

  # --- MCMC controls ---
  if (niter <= 0 || burnin < 0 || thin <= 0)
    stop("'niter' and 'thin' must be > 0, and 'burnin' >= 0.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1)
    stop("'verbose' must be a logical value (TRUE or FALSE).", call. = FALSE)
  print_progress <- if (verbose) 1L else 0L  # habilita barra/porcentaje en C++

  # --- Model frame ---
  if (is.null(data)) data <- environment(formula)
  mf <- model.frame(formula, data, na.action = NULL)
  if (anyNA(mf))
    stop("Data contains missing values; please remove or impute them.", call. = FALSE)

  y  <- model.response(mf, "numeric")
  X  <- model.matrix(attr(mf, "terms"), mf)
  coef_names <- colnames(X)
  mt <- attr(mf, "terms")

  # --- Weights ---
  w <- if (is.null(weights)) {
    rep(1, length(y))
  } else if (is.numeric(weights)) {
    if (length(weights) != length(y))
      stop("Length of 'weights' != length of response.", call. = FALSE)
    weights
  } else if (inherits(weights, "formula")) {
    model.frame(weights, data)[[1L]]
  } else stop("'weights' must be numeric or a one-sided formula.", call. = FALSE)

  p <- ncol(X)

  # --- Prior (unified -> bqr_prior interno) ---
  pri <- if (is.null(prior)) {
    as_bqr_prior(prior(), p = p, names_x = coef_names)
  } else if (inherits(prior, "prior")) {
    as_bqr_prior(prior, p = p, names_x = coef_names)
  } else if (inherits(prior, "bqr_prior")) { # compatibilidad hacia atrás
    prior
  } else {
    stop("'prior' must be NULL, a 'prior' object (see prior()), or a 'bqr_prior' (legacy).", call. = FALSE)
  }

  # --- Detectar si el usuario definió explícitamente la prior de sigma (NO confundir con defaults)
  user_defined_sigma_prior <- .user_defined_sigma(prior)

  ## --- Avisos sobre prior de sigma SOLO cuando corresponde
  if (method %in% c("score", "approximate") && isTRUE(user_defined_sigma_prior)) {
    warning(
      sprintf("Method '%s' does not estimate sigma; 'sigma_shape' and 'sigma_rate' in prior will be ignored.", method),
      call. = FALSE
    )
  }
  if (method == "ald" && isFALSE(estimate_sigma) && isTRUE(user_defined_sigma_prior)) {
    warning("With method='ald' and estimate_sigma=FALSE, 'sigma_shape' and 'sigma_rate' in prior will be ignored.",
            call. = FALSE)
  }

  # --- Pesos normalizados según método ---
  w_norm <- w / mean(w)

  # --- Chequeo de soporte 'fix_sigma' en backend ALD ---
  supports_fix_sigma <- FALSE
  if (method == "ald") {
    supports_fix_sigma <- tryCatch({
      "fix_sigma" %in% names(formals(.MCMC_BWQR_AL))
    }, error = function(e) FALSE)
  }

  # --- Backend para un tau ---
  run_backend_one <- function(tau_i) {
    draws_i <- switch(method,
                      "ald" = {
                        if (isTRUE(estimate_sigma)) {
                          .MCMC_BWQR_AL(
                            y, X, w_norm,
                            tau            = tau_i,
                            n_mcmc         = niter,
                            burnin         = burnin,
                            thin           = thin,
                            b_prior_mean   = pri$b0,
                            B_prior_prec   = solve(pri$B0),
                            c0             = pri$c0 %||% 0.001,
                            C0             = pri$C0 %||% 0.001,
                            print_progress = print_progress
                          )
                        } else {
                          if (!supports_fix_sigma) {
                            stop("Backend '.MCMC_BWQR_AL' does not support 'fix_sigma'. ",
                                 "Use estimate_sigma=TRUE or update the backend to accept fix_sigma.", call. = FALSE)
                          }
                          .MCMC_BWQR_AL(
                            y, X, w_norm,
                            tau            = tau_i,
                            n_mcmc         = niter,
                            burnin         = burnin,
                            thin           = thin,
                            b_prior_mean   = pri$b0,
                            B_prior_prec   = solve(pri$B0),
                            fix_sigma      = 1,  # σ^2 fija = 1
                            print_progress = print_progress
                          )
                        }
                      },
                      "score" = .MCMC_BWQR_SL(
                        y, X, w_norm,
                        tau            = tau_i,
                        n_mcmc         = niter,
                        burnin         = burnin,
                        thin           = thin,
                        b_prior_mean   = pri$b0,
                        B_prior_prec   = solve(pri$B0),
                        print_progress = print_progress
                      ),
                      "approximate" = .MCMC_BWQR_AP(
                        y, X, w,
                        n_mcmc         = niter,
                        burnin         = burnin,
                        thin           = thin,
                        tau            = tau_i,
                        b_prior_mean   = pri$b0,
                        B_prior_prec   = solve(pri$B0),
                        print_progress = print_progress
                      )
    )

    # --- Limpieza y coerción de draws ---
    if (is.list(draws_i) && !is.data.frame(draws_i)) {
      draws_i <- lapply(draws_i, function(x) if (is.numeric(x) || is.matrix(x)) x)
      draws_i <- draws_i[!vapply(draws_i, is.null, logical(1))]
      draws_i <- do.call(cbind, draws_i)
    }
    if (is.data.frame(draws_i))
      draws_i <- data.matrix(draws_i[vapply(draws_i, is.numeric, logical(1))])

    draws_i <- as.matrix(draws_i)
    storage.mode(draws_i) <- "numeric"
    if (anyNA(draws_i))
      stop("Backend returned non-numeric values; cannot summarize.", call. = FALSE)

    # Manejo robusto de columnas diagnósticas
    diag_cols <- c("accept_rate", "n_mcmc", "burnin", "thin", "n_samples")
    cn <- colnames(draws_i)
    keep_idx <- if (is.null(cn)) rep(TRUE, ncol(draws_i)) else !(cn %in% diag_cols)
    accept_rate_i <- if (!is.null(cn) && "accept_rate" %in% cn) mean(draws_i[, "accept_rate"]) else NA_real_
    draws_i <- draws_i[, keep_idx, drop = FALSE]

    # Nombres de columnas y betas
    if (ncol(draws_i) >= p) {
      if (is.null(colnames(draws_i))) colnames(draws_i) <- paste0("V", seq_len(ncol(draws_i)))
      colnames(draws_i)[1:p] <- coef_names
    } else if (is.null(colnames(draws_i))) {
      colnames(draws_i) <- paste0("V", seq_len(ncol(draws_i)))
    }

    beta_hat_i <- if (ncol(draws_i) >= p) colMeans(draws_i[, seq_len(p), drop = FALSE]) else numeric(0)
    names(beta_hat_i) <- if (length(beta_hat_i) > 0) coef_names else character(0)

    list(draws = draws_i, beta = beta_hat_i, accept_rate = accept_rate_i)
  }

  # --- Ejecutar para todos los taus ---
  fits <- lapply(taus, run_backend_one)
  names(fits) <- paste0("tau=", formatC(taus, format = "f", digits = 3))

  runtime <- proc.time()[["elapsed"]] - tic

  # --- Salida ---
  if (length(taus) == 1L) {
    out <- list(
      beta         = fits[[1]]$beta,
      draws        = fits[[1]]$draws,
      accept_rate  = fits[[1]]$accept_rate,
      warmup       = burnin,
      thin         = thin,
      runtime      = runtime,
      method       = method,
      quantile     = taus,
      prior        = pri,
      terms        = mt,
      model        = mf,
      formula      = formula
    )
    out$call$formula <- formula
    class(out) <- c("bwqr_fit", "bqr.svy")
    return(out)
  } else {
    beta_mat  <- do.call(cbind, lapply(fits, `[[`, "beta"))     # p × K
    colnames(beta_mat) <- names(fits)
    draws_list <- lapply(fits, `[[`, "draws")                   # list of matrices
    acc_vec    <- vapply(fits, `[[`, numeric(1), "accept_rate")
    names(acc_vec) <- names(fits)

    out <- list(
      beta           = beta_mat,
      draws          = draws_list,
      accept_rate    = acc_vec,
      warmup         = burnin,
      thin           = thin,
      runtime        = runtime,
      method         = method,
      quantile       = taus,
      prior          = pri,
      terms          = mt,
      model          = mf,
      formula        = formula,
      estimate_sigma = estimate_sigma
    )
    class(out) <- c("bwqr_fit_multi", "bqr.svy")
    return(out)
  }
}
