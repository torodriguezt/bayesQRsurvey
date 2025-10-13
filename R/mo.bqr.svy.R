# =====================================================
#  Helpers
# =====================================================

if (!exists("%||%"))
  `%||%` <- function(a, b) if (is.null(a)) b else a

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

  if (inherits(pr, "mo_bqr_prior")) {
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

# =====================================================
#  MODEL FITTER (multiple-output/multivariate)
# =====================================================

#' Multiple-Output Bayesian quantile regression for complex survey data
#'
#' mo.bqr.svy implements a Bayesian approach to multiple-output quantile regression
#' for complex survey data analysis. The method builds a quantile region based on
#' a directional approach. To improve computational efficiency, an Expectation-Maximization (EM)
#' algorithm is implemented instead of the usual Markov Chain Monte Carlo (MCMC).
#'
#' @param formula a symbolic description of the model to be fit.
#' @param weights an optional numerical vector containing the survey weights. If \code{NULL}, equal weights are used.
#' @param data an optional data frame containing the variables in the model.
#' @param quantile numerical scalar or vector containing quantile(s) of interest (default=0.5).
#' @param prior a \code{bqr_prior} object of class "prior". If omitted, a vague prior is assumed (see \code{\link{prior}}).
#' @param U an optional \eqn{d \times K}-matrix of directions, where \eqn{d} indicates the response variable dimension
#' and \eqn{K} indicates indicates the number of directions.
#' @param gamma_U an optional list with length equal to \eqn{K} for which each element corresponds to
#' \eqn{d \times (d-1)}-matrix of ortoghonal basis for each row of \code{U}.
#' @param n_dir numerical scalar corresponding to the number of directions (if \code{U} and \code{gamma_U} are not supplied).
#' @param epsilon numerical scalar indicating the convergence tolerance for the EM algorithm (default = 1e-6).
#' @param max_iter numerical scalar indicating maximum number of EM iterations (default = 1000).
#' @param verbose logical flag indicating whether to print progress messages (default=FALSE).
#' @param estimate_sigma logical flag indicating whether to estimate the scale parameter
#' when method = "ald" (default=FALSE and \eqn{\sigma^2} is set to 1)
#'
#' @return An object of class \code{"mo.bqr.svy"} containing:
#'   \item{call}{The matched call}
#'   \item{formula}{The model formula}
#'   \item{terms}{The terms object}
#'   \item{quantile}{Vector of fitted quantiles}
#'   \item{prior}{List of priors used for each quantile}
#'   \item{fit}{List of fitted results for each quantile, each containing one sub-list per direction}
#'   \item{coefficients}{Coefficients organized by quantile}
#'   \item{sigma}{List of scale parameters by quantile and direction.
#'                If \code{estimate_sigma = FALSE}, all entries are fixed at 1.
#'                If \code{estimate_sigma = TRUE}, each entry contains the
#'                estimated value of \eqn{\sigma} (posterior mode from EM).}
#'   \item{n_dir}{Number of directions}
#'   \item{U}{Matrix of projection directions (\eqn{d \times K})}
#'   \item{Gamma_list}{List of orthogonal complement bases, one per direction}
#'   \item{n_obs}{Number of observations}
#'   \item{n_vars}{Number of covariates}
#'   \item{response_dim}{Dimension of the response \eqn{d}}
#'   \item{estimate_sigma}{Logical flag indicating whether the scale parameter
#'                         \eqn{\sigma^2} was estimated (\code{TRUE}) or fixed at 1 (\code{FALSE}).}
#'
#' @references
#' Nascimento, M. L. & \enc{Gonçalves}{Goncalves}, K. C. M. (2024).
#' Bayesian Quantile Regression Models for Complex Survey Data Under Informative Sampling.
#' \emph{Journal of Survey Statistics and Methodology}, 12(4), 1105–1130.
#' <doi:10.1093/jssam/smae015>
#'
#' @examples
#' library(MASS)
#'
#' # Generate population data
#' set.seed(123)
#' N    <- 10000
#' data <- mvrnorm(N,rep(0,3),matrix(c(4,0,2,0,1,1.5,2,1.5,9),3,3))
#' x_p  <- as.matrix(data[,1])
#' y_p  <- data[,2:3]+cbind(rep(0,N),x_p)
#'
#' # Generate sample data
#' n <- 500
#' z_aux <- rnorm(N, mean = 1 + y_p, sd=.5)
#' p_aux <- 1 / (1 + exp(2.5 - 0.5 * z_aux))
#' s_ind <- sample(1:N, n, replace = FALSE, prob = p_aux)
#' y_s <- y_p[s_ind,]
#' x_s <- x_p[s_ind,]
#' w   <- 1 / p_aux[s_ind]
#' data_s<- data.frame(y1 = y_s[,1], y2 = y_s[,2],  x1 = x_s, w = w)
#'
#' # Basic usage with default priors when U and gamma_U are given
#' fit1 <- mo.bqr.svy(cbind(y1, y2) ~ x1, weights = w, data = data_s, quantile = c(0.1, 0.2),
#'  U = matrix(c(0,1,1/sqrt(2),1/sqrt(2)),2),
#'  gamma_U = list(c(1,0),c(1/sqrt(2),-1/sqrt(2))))
#'
#' # Basic usage with default priors when n_dir is given
#' fit2 <- mo.bqr.svy(cbind(y1, y2) ~ x1, weights = w, data = data_s, quantile = c(0.1, 0.2), n_dir=2)
#'
#' @export
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom pracma nullspace
mo.bqr.svy <- function(formula,
                       weights  = NULL,
                       data     = NULL,
                       quantile = 0.5,
                       prior    = NULL,
                       U        = NULL,
                       gamma_U  = NULL,
                       n_dir    = NULL,
                       epsilon  = 1e-6,
                       max_iter = 1000,
                       verbose  = FALSE,
                       estimate_sigma = FALSE) {

  algorithm <- "em"
  if (missing(data)) stop("'data' must be provided.")

  # --- quantiles ---
  if (length(quantile) == 0) stop("'quantile' cannot be empty.")
  if (any(!is.finite(quantile))) stop("'quantile' must be numeric and finite.")
  if (any(quantile <= 0 | quantile >= 1)) stop("'quantile' must be in (0,1).")
  quantile <- sort(unique(quantile))

  # --- model frame ---
  mf <- model.frame(formula, data)
  y  <- model.response(mf)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  if (!is.matrix(y)) stop("'y' must be a numeric matrix or vector.")
  if (any(!is.finite(y))) stop("Response 'y' contains non-finite values.")
  n <- nrow(y); d <- ncol(y)

  X  <- model.matrix(attr(mf, "terms"), mf)
  if (nrow(X) != n) stop("nrow(X) must match nrow(y).")
  p  <- ncol(X)
  coef_names <- colnames(X)

  # --- weights ---
  wts <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  if (length(wts) != n)  stop("'weights' must have length n.")
  if (any(!is.finite(wts)) || any(wts <= 0)) stop("Invalid weights.")
  wts <- wts / mean(wts)

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package 'pracma' is required for nullspace calculation")
  }

  # ---------- Build U and Gamma_list ----------
  if (!is.null(U)) {
    U <- as.matrix(U)
    if (nrow(U) != d) stop("U must have d rows.")
    K <- ncol(U)

    # normalize U columns
    tol_u <- 1e-12
    for (k in seq_len(K)) {
      nk <- sqrt(sum(U[, k]^2))
      if (!is.finite(nk) || nk <= tol_u)
        stop(sprintf("U[, %d] must be a non-zero finite vector.", k))
      if (abs(nk - 1) > 1e-8) {
        if (verbose) message(sprintf("Normalizing U[, %d] (||u||=%.6g).", k, nk))
        U[, k] <- U[, k] / nk
      }
    }

    if (!is.null(gamma_U)) {
      if (!is.list(gamma_U)) stop("'gamma_U' must be a list of K matrices.")
      if (length(gamma_U) != K) stop("length(gamma_U) must equal ncol(U) = K.")

      Gamma_list <- vector("list", K)
      for (k in seq_len(K)) {
        if (d == 1) {
          Gamma_list[[k]] <- matrix(numeric(0), 1, 0)
        } else {
          Gk <- as.matrix(gamma_U[[k]])
          if (nrow(Gk) != d) stop(sprintf("gamma_U[[%d]] must have %d rows.", k, d))
          if (ncol(Gk) != (d - 1))
            stop(sprintf("gamma_U[[%d]] must have %d columns (d-1).", k, d - 1))
          u_k <- U[, k]
          tol <- 1e-8
          ortho_err <- max(abs(drop(crossprod(Gk, u_k))))
          if (ortho_err > tol)
            stop(sprintf("gamma_U[[%d]] is not orthogonal to U[,%d].",
                         k, k, ortho_err, tol))
          if (qr(Gk)$rank < (d - 1))
            stop(sprintf("gamma_U[[%d]] does not have full column rank d-1.", k))
          if (verbose) {
            gram_err <- max(abs(crossprod(Gk) - diag(d - 1)))
            if (gram_err > 1e-6) {
              message(sprintf("Warning: gamma_U[[%d]] columns are not perfectly orthonormal (max Gram error %.2e). Proceeding as provided.", k, gram_err))
            }
          }
          Gamma_list[[k]] <- Gk
        }
      }

    } else {
      if (d == 1) {
        Gamma_list <- replicate(K, matrix(numeric(0), 1, 0), simplify = FALSE)
      } else {
        Gamma_list <- lapply(seq_len(K), function(k) pracma::nullspace(t(U[, k])))
      }
    }

  } else {
    if (!is.null(gamma_U)) stop("If you provide 'gamma_U', you must also provide 'U'.")
    if (d == 1) {
      U <- matrix(1.0, 1, 1); K <- 1
      Gamma_list <- list(matrix(numeric(0), 1, 0))
    } else {
      K <- max(1L, as.integer(n_dir))
      U <- matrix(NA_real_, d, K)
      Gamma_list <- vector("list", K)
      for (k in seq_len(K)) {
        u_k <- rnorm(d); u_k <- u_k / sqrt(sum(u_k^2))
        U[, k] <- u_k
        Gamma_list[[k]] <- pracma::nullspace(t(u_k))
      }
    }
  }

  # ---------- Unified prior ----------
  pri <- if (is.null(prior)) {
    as_mo_bqr_prior(prior(), p = p, d = d, names_x = coef_names, names_y = NULL)
  } else if (inherits(prior, "prior")) {
    as_mo_bqr_prior(prior, p = p, d = d, names_x = coef_names, names_y = NULL)
  } else if (inherits(prior, "mo_bqr_prior")) {
    prior
  } else {
    stop("'prior' must be NULL, a 'prior' object (see prior()), or a 'mo_bqr_prior' (legacy).", call. = FALSE)
  }

  user_defined_sigma_prior <- .user_defined_sigma(prior)
  if (isFALSE(estimate_sigma) && isTRUE(user_defined_sigma_prior)) {
    warning("With estimate_sigma=FALSE, 'sigma_shape' and 'sigma_rate' in prior will be ignored.", call. = FALSE)
  }

  # --- Backend 'fix_sigma' support?
  backend_supports_fix <- tryCatch({
    "fix_sigma" %in% names(formals(.bwqr_weighted_em_cpp_sep))
  }, error = function(e) FALSE)

  results <- vector("list", length(quantile))
  names(results) <- paste0("tau=", formatC(quantile, digits = 3, format = "f"))

  for (qi in seq_along(quantile)) {
    qtau <- quantile[qi]
    pr   <- pri

    if (verbose) message(sprintf("Fitting tau = %.3f (d=%d, K=%d)", qtau, d, K))

    direction_results <- vector("list", K)
    names(direction_results) <- paste0("dir_", seq_len(K))

    for (k in seq_len(K)) {
      u_k <- U[, k]
      gamma_uk <- Gamma_list[[k]]
      r_k <- if (d > 1) ncol(gamma_uk) else 0
      p_ext <- p + r_k

      # prior blocks
      beta_mean_k <- pr$beta_mean
      beta_cov_k  <- pr$beta_cov

      if (r_k > 0) {
        if (is.null(pr$beta_star_mean)) {
          gamma_mean_k <- rep(0, r_k)
        } else if (length(pr$beta_star_mean) == 1L) {
          gamma_mean_k <- rep(pr$beta_star_mean, r_k)
        } else if (length(pr$beta_star_mean) == r_k) {
          gamma_mean_k <- pr$beta_star_mean
        } else {
          stop(sprintf("Length of prior$beta_star_mean (%d) must be 1 or r_k=%d.",
                       length(pr$beta_star_mean), r_k))
        }

        if (is.null(pr$beta_star_cov)) {
          gamma_cov_k <- diag(1e6, r_k)
        } else if (is.numeric(pr$beta_star_cov) && length(pr$beta_star_cov) == 1L) {
          gamma_cov_k <- diag(as.numeric(pr$beta_star_cov), r_k)
        } else if (is.numeric(pr$beta_star_cov) && is.null(dim(pr$beta_star_cov)) &&
                   length(pr$beta_star_cov) == r_k) {
          gamma_cov_k <- diag(as.numeric(pr$beta_star_cov), r_k)
        } else if (is.matrix(pr$beta_star_cov) && all(dim(pr$beta_star_cov) == c(r_k, r_k))) {
          gamma_cov_k <- pr$beta_star_cov
        } else {
          stop(sprintf("prior$beta_star_cov must be scalar, length r_k vector, or r_k x r_k matrix (here r_k=%d).", r_k))
        }
      } else {
        gamma_mean_k <- numeric(0)
        gamma_cov_k  <- matrix(numeric(0), 0, 0)
      }

      mu0_ext <- c(beta_mean_k, gamma_mean_k)
      sigma0_ext <- matrix(0, p_ext, p_ext)
      sigma0_ext[1:p, 1:p] <- beta_cov_k
      if (r_k > 0) sigma0_ext[(p+1):p_ext, (p+1):p_ext] <- gamma_cov_k

      u_k_matrix     <- matrix(u_k, ncol = 1)
      gamma_k_matrix <- if (d > 1) gamma_uk else matrix(numeric(0), d, 0)

      # --- Backend call: estimate sigma or fix to 1
      cpp_result <- if (isTRUE(estimate_sigma)) {
        .bwqr_weighted_em_cpp_sep(
          y        = y,
          x        = X,
          w        = wts,
          u        = u_k_matrix,
          gamma_u  = gamma_k_matrix,
          tau      = qtau,
          mu0      = mu0_ext,
          sigma0   = sigma0_ext,
          a0       = pr$sigma_shape,
          b0       = pr$sigma_rate,
          eps      = epsilon,
          max_iter = max_iter,
          verbose  = verbose
        )
      } else if (backend_supports_fix) {
        .bwqr_weighted_em_cpp_sep(
          y        = y,
          x        = X,
          w        = wts,
          u        = u_k_matrix,
          gamma_u  = gamma_k_matrix,
          tau      = qtau,
          mu0      = mu0_ext,
          sigma0   = sigma0_ext,
          a0       = pr$sigma_shape,
          b0       = pr$sigma_rate,
          eps      = epsilon,
          max_iter = max_iter,
          verbose  = verbose,
          fix_sigma = 1.0
        )
      } else {
        a0_use <- 1e9; b0_use <- a0_use + 1
        if (verbose) message("Backend without 'fix_sigma'; using highly concentrated prior at sigma = 1.")
        out <- .bwqr_weighted_em_cpp_sep(
          y        = y,
          x        = X,
          w        = wts,
          u        = u_k_matrix,
          gamma_u  = gamma_k_matrix,
          tau      = qtau,
          mu0      = mu0_ext,
          sigma0   = sigma0_ext,
          a0       = a0_use,
          b0       = b0_use,
          eps      = epsilon,
          max_iter = max_iter,
          verbose  = verbose
        )
        if (!is.null(out$sigma)) out$sigma[] <- 1.0
        out
      }

      beta_k  <- as.numeric(cpp_result$beta[1, ])
      sigma_k <- if (isTRUE(estimate_sigma)) as.numeric(cpp_result$sigma[1]) else 1.0
      covariate_names <- c(coef_names,
                           if (r_k > 0) paste0("gamma_", seq_len(r_k)) else character(0))
      names(beta_k) <- covariate_names

      direction_results[[k]] <- list(
        beta        = beta_k,
        sigma       = sigma_k,
        iter        = cpp_result$iter,
        converged   = cpp_result$converged,
        u           = u_k,
        gamma_u     = gamma_uk
      )
    }

    results[[qi]] <- list(
      directions  = direction_results,
      prior       = pri,
      U           = U,
      Gamma_list  = Gamma_list,
      quantile    = qtau
    )
  }

  # --- Coefficients: list per tau (matrix coef x directions) ---
  coefficients_list <- vector("list", length(quantile))
  names(coefficients_list) <- names(results)
  for (qi in seq_along(results)) {
    dir_list <- results[[qi]]$directions
    all_rows <- unique(unlist(lapply(dir_list, function(dd) names(dd$beta))))
    mat <- matrix(NA_real_, nrow = length(all_rows), ncol = length(dir_list))
    rownames(mat) <- all_rows
    colnames(mat) <- paste0("dir", seq_along(dir_list))
    for (k in seq_along(dir_list)) {
      b <- dir_list[[k]]$beta
      mat[names(b), k] <- b
    }
    coefficients_list[[qi]] <- mat
  }

  # --- Sigma: list per tau (numeric vector per directions) ---
  sigma_list <- vector("list", length(quantile))
  names(sigma_list) <- names(results)
  for (qi in seq_along(results)) {
    dir_list <- results[[qi]]$directions
    vals <- vapply(dir_list, function(dd) {
      v <- tryCatch(as.numeric(dd$sigma)[1], error = function(e) NA_real_)
      ifelse(is.finite(v), v, NA_real_)
    }, numeric(1))
    names(vals) <- paste0("dir", seq_along(dir_list))
    sigma_list[[qi]] <- vals
  }

  structure(list(
    call            = match.call(),
    formula         = formula,
    terms           = attr(mf, "terms"),
    quantile        = quantile,
    prior           = pri,
    fit             = results,
    coefficients    = coefficients_list,
    sigma           = sigma_list,
    n_dir           = K,
    U               = U,
    Gamma_list      = Gamma_list,
    n_obs           = n,
    n_vars          = p,
    response_dim    = d,
    estimate_sigma  = isTRUE(estimate_sigma)
  ), class = "mo.bqr.svy")
}
