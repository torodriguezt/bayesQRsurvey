# ======================================================================
# Utilities & imports (posterior-based summaries/diagnostics)
# ======================================================================

#' @keywords internal
#' @importFrom posterior summarize_draws quantile2
NULL

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.atomic(a) && length(a) == 1L) return(if (is.na(a)) b else a)
  a
}

# Internal: posterior summaries with consistent columns
summarise_draws_custom <- function(draws, probs = c(0.025, 0.975), ...) {
  if (!is.matrix(draws)) draws <- as.matrix(draws)
  if (is.null(colnames(draws))) colnames(draws) <- paste0("V", seq_len(ncol(draws)))

  s <- posterior::summarize_draws(
    draws,
    "mean", "median", "sd",
    "rhat", "ess_bulk", "ess_tail",
    ~posterior::quantile2(.x, probs = probs)
  )

  # Add lower/upper_ci when exactly two probs are provided
  if (length(probs) == 2) {
    if (all(c("q2.5", "q97.5") %in% names(s))) {
      s$lower_ci <- s$`q2.5`
      s$upper_ci <- s$`q97.5`
    } else {
      qcols <- grep("^q[0-9]", names(s), value = TRUE)
      if (length(qcols) >= 2) {
        qnum <- suppressWarnings(as.numeric(sub("^q", "", sub("\\.#", ".", qcols))))
        ord <- order(qnum)
        s$lower_ci <- s[[qcols[ord[1]]]]
        s$upper_ci <- s[[qcols[ord[length(ord)]]]]
      }
    }
  }

  wanted <- c("variable","mean","median","sd","rhat","ess_bulk","ess_tail",
              "q2.5","q97.5","lower_ci","upper_ci")
  s[, intersect(wanted, names(s)), drop = FALSE]
}


# ======================================================================
# ANCHOR: SUMMARY (single help page "summary.bayesQRsurvey")
# ======================================================================

#' Summary methods for bayesQRsurvey
#'@description
#' summary.bayesQRsurvey is an S3 method that summarizes the output of the
#' \code{bqr.svy} or \code{mo.bqr.svy} function. For the \code{bqr.svy} the posterior mean,
#' posterior credible interval and convergence diagnostics are calculated. For the \code{mo.bqr.svy}
#' the iterations for convergence, the MAP and the direction are calculated.
#' @name summary.bayesQRsurvey
#' @docType methods
NULL


#' @rdname summary.bayesQRsurvey
#' @title Summary of \code{bqr.svy} fits
#' @param object An object of class \code{bqr.svy}.
#' @param probs Two-element numeric vector with credible interval probabilities.
#'   Default \code{c(0.025, 0.975)}.
#' @param digits Integer; number of decimals used by printing helpers.
#' @param ... Unused.
#' @return An object of class \code{summary.bqr.svy} with one block per \eqn{\tau}.
#' @exportS3Method summary bqr.svy
summary.bqr.svy <- function(object, probs = c(0.025, 0.975), digits = 3, ...) {
  stopifnot(inherits(object, "bqr.svy"))

  meta <- list(
    n_chains    = object$n_chains %||% 1L,
    warmup      = object$warmup   %||% 0L,
    thin        = object$thin     %||% 1L,
    accept_rate = object$accept_rate %||% NA_real_,
    runtime     = object$runtime  %||% NA_real_
  )

  ic_by_name <- function(D, vars, probs) {
    lower <- upper <- rep(NA_real_, length(vars))
    cn <- colnames(D)
    for (i in seq_along(vars)) {
      j <- match(vars[i], cn)
      if (!is.na(j)) {
        vcol <- D[, j]
        lower[i] <- stats::quantile(vcol, probs = probs[1], na.rm = TRUE)
        upper[i] <- stats::quantile(vcol, probs = probs[2], na.rm = TRUE)
      }
    }
    list(lower = lower, upper = upper)
  }

  make_block <- function(D, tau) {
    D <- if (is.data.frame(D)) data.matrix(D) else as.matrix(D)
    storage.mode(D) <- "numeric"
    if (!nrow(D) || !ncol(D)) stop("Empty draws matrix.", call. = FALSE)

    stats <- summarise_draws_custom(D, probs = probs)

    # Include sigma in main table only when it was estimated
    if (isTRUE(object$estimate_sigma)) {
      coef_rows <- stats
    } else {
      coef_rows <- stats[stats$variable != "sigma", , drop = FALSE]
    }

    # Main table: estimate + credible interval (3 numeric columns)
    main_cols    <- c("variable", "mean", "lower_ci", "upper_ci")
    coef_summary <- coef_rows[, intersect(main_cols, names(coef_rows)), drop = FALSE]

    # Re-compute CI directly from draws for accuracy
    vars <- coef_summary$variable
    ic   <- ic_by_name(D, vars, probs)
    coef_summary$lower_ci <- ic$lower
    coef_summary$upper_ci <- ic$upper

    list(
      tau          = tau,
      coef_summary = coef_summary,
      full_summary = stats,
      n_draws      = nrow(D),
      meta         = meta,
      probs        = probs,
      digits       = digits
    )
  }

  if (is.list(object$draws)) {
    per_tau <- Map(make_block, object$draws, object$quantile)
    names(per_tau) <- paste0("tau=", formatC(object$quantile, format = "f", digits = 3))
  } else {
    per_tau <- list(make_block(object$draws, object$quantile))
    names(per_tau) <- paste0("tau=", formatC(object$quantile, format = "f", digits = 3))
  }

  res <- list(
    call      = object$call %||% NULL,
    method    = object$method %||% object$algorithm %||% NA_character_,
    quantiles = object$quantile,
    per_tau   = per_tau
  )
  class(res) <- "summary.bqr.svy"
  res
}


#' @rdname summary.bayesQRsurvey
#' @title Summary of \code{mo.bqr.svy} fits
#' @param object An object of class \code{mo.bqr.svy}.
#' @param digits Integer; number of decimals used by printing helpers. Default \code{4}.
#' @param ... Unused.
#' @return An object of class \code{summary.mo.bqr.svy}.
#' @exportS3Method summary mo.bqr.svy
summary.mo.bqr.svy <- function(object, digits = 4, ...) {
  stopifnot(inherits(object, "mo.bqr.svy"))
  `%||%` <- function(a, b) if (is.null(a)) b else a

  one_dir <- function(dir_res, tau, k) {
    beta <- dir_res$beta
    coef_tab <- data.frame(
      parameter   = names(beta),
      MAP         = as.numeric(beta),
      check.names = FALSE
    )
    list(
      tau       = tau,
      dir_id    = k,
      coef_tab  = coef_tab,
      sigma     = as.numeric(dir_res$sigma),
      iter      = dir_res$iter %||% NA_integer_,
      converged = isTRUE(dir_res$converged)
    )
  }

  blocks <- list()
  for (qi in seq_along(object$quantile)) {
    tau  <- object$quantile[qi]
    dirs <- object$fit[[qi]]$directions
    for (k in seq_along(dirs)) {
      blocks[[length(blocks) + 1]] <- one_dir(dirs[[k]], tau, k)
    }
  }

  out <- list(
    call         = object$call %||% NULL,
    algorithm    = object$algorithm %||% "em",
    tolerance    = object$tolerance %||% object$epsilon %||% NA_real_,
    quantiles    = object$quantile,
    n_dir        = object$n_dir,
    n_obs        = object$n_obs,
    n_vars       = object$n_vars,
    response_dim = object$response_dim,
    digits       = digits,
    blocks       = blocks
  )
  class(out) <- "summary.mo.bqr.svy"
  out
}


# (Hidden) print methods for summary objects: exported, no Rd
#' @noRd
#' @exportS3Method print summary.bqr.svy
print.summary.bqr.svy <- function(x, ...) {
  stopifnot(inherits(x, "summary.bqr.svy"))
  cat("\nMethod: ", x$method %||% "NA", "\n", sep = "")
  cat("Quantiles: ",
      paste(formatC(x$quantiles, format = "f", digits = 3), collapse = ", "),
      "\n\n", sep = "")

  acc_global <- tryCatch({
    vals <- vapply(x$per_tau, function(b) {
      a <- b$meta$accept_rate
      if (length(a) == 0) NA_real_ else mean(a, na.rm = TRUE)
    }, numeric(1))
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, error = function(e) NA_real_)

  for (nm in names(x$per_tau)) {
    blk <- x$per_tau[[nm]]
    cat("== ", nm, " ==\n", sep = "")

    acc_here <- blk$meta$accept_rate
    acc_here <- if (length(acc_here) <= 1L) as.numeric(acc_here) else
      suppressWarnings(mean(acc_here, na.rm = TRUE))

    if (is.finite(acc_here)) {
      cat("  Acceptance rate (avg): ", round(acc_here, 3), "\n", sep = "")
    } else if (is.finite(acc_global)) {
      cat("  Acceptance rate (avg): ", round(acc_global, 3), "\n", sep = "")
    }

    cat("  Draws: ", blk$n_draws,
        " | Warmup: ", blk$meta$warmup,
        " | Thin: ", blk$meta$thin, "\n\n", sep = "")

    cs <- blk$coef_summary
    show_cols <- c("variable", "mean", "lower_ci", "upper_ci")
    show_cols <- intersect(show_cols, colnames(cs))
    if (!length(show_cols)) {
      print(cs)
    } else {
      df <- cs[, show_cols, drop = FALSE]
      num_cols <- setdiff(show_cols, "variable")
      for (cl in num_cols) {
        df[[cl]] <- ifelse(is.finite(df[[cl]]), round(df[[cl]], blk$digits), df[[cl]])
      }
      print(df, row.names = FALSE)
    }
    cat("\n")
  }
  invisible(x)
}

#' @noRd
#' @exportS3Method print summary.mo.bqr.svy
print.summary.mo.bqr.svy <- function(x, max_dir = 8, ...) {
  stopifnot(inherits(x, "summary.mo.bqr.svy"))
  dig  <- x$digits
  rule <- strrep("\u2500", 56)

  # --- Header ---
  cat("\n")
  cat("  Multiple-Output Bayesian Quantile Regression (Summary)\n")
  cat("  ", rule, "\n", sep = "")

  qs <- paste(formatC(x$quantiles, format = "f", digits = 3), collapse = ", ")
  cat("  Quantiles  : ", qs, "\n", sep = "")
  cat("  Directions : ", x$n_dir, "\n", sep = "")
  if (!is.null(x$n_obs))
    cat("  Sample     : ", x$n_obs, " obs", sep = "")
  if (!is.null(x$response_dim))
    cat(", ", x$response_dim, " responses", sep = "")
  cat("\n")
  cat("  Estimation : EM (posterior mode / MAP)\n")
  cat("  ", rule, "\n", sep = "")

  # --- Order blocks by tau, then direction ---
  ord <- order(
    vapply(x$blocks, function(b) b$tau,    numeric(1)),
    vapply(x$blocks, function(b) b$dir_id, numeric(1))
  )
  blocks <- x$blocks[ord]

  # --- Group blocks by tau ---
  tau_vals <- unique(vapply(blocks, function(b) b$tau, numeric(1)))

  for (tau in tau_vals) {
    tau_blocks <- blocks[vapply(blocks, function(b) b$tau == tau, logical(1))]
    n_dir <- length(tau_blocks)
    tau_lab <- formatC(tau, format = "f", digits = 3)

    # Build coefficient matrix (rows = params, cols = dirs)
    param_names <- tau_blocks[[1]]$coef_tab$parameter
    coef_mat <- matrix(NA_real_, nrow = length(param_names), ncol = n_dir)
    rownames(coef_mat) <- param_names
    colnames(coef_mat) <- paste0("dir", vapply(tau_blocks, function(b) b$dir_id, numeric(1)))
    for (j in seq_along(tau_blocks))
      coef_mat[, j] <- tau_blocks[[j]]$coef_tab$MAP

    cat("\n  Coefficients (tau = ", tau_lab, ")\n", sep = "")
    .print_coef_table(coef_mat, digits = dig, max_dir = max_dir, indent = 4)

    # --- Convergence summary ---
    iters  <- vapply(tau_blocks, function(b)
      if (!is.na(b$iter)) b$iter else NA_integer_, numeric(1))
    conv   <- vapply(tau_blocks, function(b) isTRUE(b$converged), logical(1))
    sigmas <- vapply(tau_blocks, function(b) b$sigma, numeric(1))

    n_conv <- sum(conv)
    n_tot  <- length(conv)
    valid_iters <- iters[!is.na(iters)]

    cat("\n  EM convergence: ", n_conv, "/", n_tot, " directions converged\n", sep = "")
    if (length(valid_iters) > 0)
      cat("  Iterations    : min = ", min(valid_iters),
          ", median = ", stats::median(valid_iters),
          ", max = ", max(valid_iters), "\n", sep = "")

    # Show non-converged directions if any
    if (n_conv < n_tot) {
      not_conv <- which(!conv)
      cat("  Not converged : dir ",
          paste(not_conv, collapse = ", "), "\n", sep = "")
    }

    # Sigma summary
    unique_sig <- unique(round(sigmas, dig))
    if (length(unique_sig) == 1L) {
      cat("  Sigma         : ", unique_sig, " (all directions)\n", sep = "")
    } else {
      cat("  Sigma         : min = ", min(sigmas),
          ", max = ", max(sigmas), "\n", sep = "")
    }
  }

  cat("\n")
  invisible(x)
}


# Internal helper: tidy any summary object into a common schema (no bwqr)
..tidy_bayesQRsurvey_summary <- function(s, target_tau = 0.5) {
  if (inherits(s, "summary.bqr.svy")) {
    taus <- vapply(s$per_tau, `[[`, numeric(1), "tau")
    idx  <- which.min(abs(taus - target_tau))
    df   <- s$per_tau[[idx]]$coef_summary
  } else if (inherits(s, "summary.mo.bqr.svy")) {
    blocks <- s$blocks
    make_df <- function(b) {
      data.frame(
        variable = paste0(
          b$coef_tab$parameter,
          " [dir ", b$dir_id, ", tau=", formatC(b$tau, format = "f", digits = 3), "]"
        ),
        mean     = b$coef_tab$MAP,
        sd       = NA_real_,
        rhat     = NA_real_,
        ess_bulk = NA_real_,
        ess_tail = NA_real_,
        lower_ci = NA_real_,
        upper_ci = NA_real_,
        check.names = FALSE
      )
    }
    df <- do.call(rbind, lapply(blocks, make_df))
  } else {
    stop("Unknown summary object class: ", paste(class(s), collapse = ", "))
  }

  wanted <- c("variable","mean","sd","rhat","ess_bulk","ess_tail","lower_ci","upper_ci")
  miss <- setdiff(wanted, names(df))
  for (nm in miss) df[[nm]] <- NA_real_
  df[, wanted, drop = FALSE]
}

# (Hidden) combine multiple fits into a comparative table
#' @noRd
#' @exportS3Method summary list
summary.list <- function(object, ..., methods = NULL, target_tau = 0.5, digits = 3) {
  supported <- c("bqr.svy","mo.bqr.svy")
  is_supported <- function(x) any(inherits(x, supported))
  if (!length(object) || !all(vapply(object, is_supported, logical(1)))) {
    return(NextMethod())
  }

  smry <- lapply(object, function(f) summary(f, ...))

  if (is.null(methods)) {
    methods <- vapply(object, function(f)
      (f$method %||% class(f)[1]) %||% "model", character(1))
  }
  if (length(methods) != length(smry)) {
    stop("'methods' must be NULL or have the same length as 'object'.")
  }

  tidied <- Map(function(s, lab) {
    df <- ..tidy_bayesQRsurvey_summary(s, target_tau = target_tau)
    df$Method <- lab
    df
  }, smry, methods)

  table <- do.call(rbind, tidied)
  table <- table[, c("Method", setdiff(names(table), "Method")), drop = FALSE]

  out <- list(
    table      = table,
    target_tau = target_tau,
    components = smry,
    digits     = digits
  )
  class(out) <- "summary.bayesQRsurvey"
  out
}


# ======================================================================
# ANCHOR: PRINT (single help page "print.bayesQRsurvey")
# ======================================================================

#' Print methods for bayesQRsurvey model objects
#'
#' @description
#' \code{print.bayesQRsurvey} is an S3 method that prints the content of an S3 object of class
#' \code{bqr.svy} or \code{mo.bqr.svy} to the console.
#'
#' @param x An object of class \code{"bqr.svy"} or \code{"mo.bqr.svy"},
#'   returned by \code{\link{bqr.svy}} or \code{\link{mo.bqr.svy}}.
#' @param digits Integer specifying the number of decimal places to print. Defaults to \code{3}.
#' @param ... Additional arguments that are passed to the generic \code{print()} function.
#'
#' @name print.bayesQRsurvey
#' @docType methods
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' N    <- 10000
#' x1_p <- runif(N, -1, 1)
#' x2_p <- runif(N, -1, 1)
#' y_p  <- 2 + 1.5 * x1_p - 0.8 * x2_p + rnorm(N)
#'
#' # Generate sample data
#' n <- 500
#' z_aux <- rnorm(N, mean = 1 + y_p, sd = .5)
#' p_aux <- 1 / (1 + exp(2.5 - 0.5 * z_aux))
#' s_ind <- sample(1:N, n, replace = FALSE, prob = p_aux)
#' y_s   <- y_p[s_ind]
#' x1_s  <- x1_p[s_ind]
#' x2_s  <- x2_p[s_ind]
#' w     <- 1 / p_aux[s_ind]
#' data  <- data.frame(y = y_s, x1 = x1_s, x2 = x2_s, w = w)
#'
#' # Fit a model
#' fit1 <- bqr.svy(y ~ x1 + x2, weights = w, data = data,
#'                 niter = 2000, burnin = 500, thin = 2)
#'
#' print(fit1)
#' }
NULL


#' @rdname print.bayesQRsurvey
#' @title Print a \code{bqr.svy} model
#' @exportS3Method print bqr.svy
print.bqr.svy <- function(x, digits = 3, ...) {
  cat("Bayesian Quantile Regression for survey data\n")
  cat("Method   :", x$method, "\n")
  cat("Quantile :", paste(formatC(x$quantile, digits = 3, format = "f"),
                          collapse = " "), "\n")
  cat("Formula  :", deparse(x$formula), "\n")
  if (!is.null(x$runtime)) cat("Runtime  :", round(x$runtime, 3), "sec\n")

  # Build coefficient display, appending sigma row when estimated
  beta_display <- x$beta
  if (identical(x$method, "ald") && isTRUE(x$estimate_sigma)) {
    sigma_mean <- function(m) {
      if (is.matrix(m) && "sigma" %in% colnames(m)) mean(m[, "sigma"], na.rm = TRUE)
      else NA_real_
    }
    sig_vec <- if (is.matrix(x$draws) && length(x$quantile) == 1L) {
      sigma_mean(x$draws)
    } else if (is.list(x$draws)) {
      vapply(x$draws, sigma_mean, numeric(1))
    } else {
      NA_real_
    }
    if (is.matrix(beta_display)) {
      sigma_row         <- matrix(sig_vec, nrow = 1, dimnames = list("sigma", colnames(beta_display)))
      beta_display      <- rbind(beta_display, sigma_row)
  # Build coefficient display, appending sigma row when estimated
  beta_display <- x$beta
  if (identical(x$method, "ald") && isTRUE(x$estimate_sigma)) {
    sigma_mean <- function(m) {
      if (is.matrix(m) && "sigma" %in% colnames(m)) mean(m[, "sigma"], na.rm = TRUE)
      else NA_real_
    }
    sig_vec <- if (is.matrix(x$draws) && length(x$quantile) == 1L) {
      sigma_mean(x$draws)
    } else if (is.list(x$draws)) {
      vapply(x$draws, sigma_mean, numeric(1))
    } else {
      NA_real_
    }
    if (is.matrix(beta_display)) {
      sigma_row         <- matrix(sig_vec, nrow = 1, dimnames = list("sigma", colnames(beta_display)))
      beta_display      <- rbind(beta_display, sigma_row)
    } else {
      beta_display <- c(beta_display, sigma = sig_vec)
    }
  } else if (identical(x$method, "ald") && !isTRUE(x$estimate_sigma)) {
    cat("  (sigma fixed at 1)\n")
  }

  cat("\nCoefficients (posterior means):\n")
  print(round(beta_display, digits))

  # ---- Acceptance rate ----
  if (!is.null(x$accept_rate)) {
    if (is.numeric(x$accept_rate) && length(x$accept_rate) == 1L) {
      cat("\nAcceptance rate:", round(x$accept_rate, 3), "\n")
    } else if (is.numeric(x$accept_rate) && length(x$accept_rate) > 1L) {
      cat("\nAcceptance rate by quantile:\n")
      print(round(x$accept_rate, 3))
    }
  }

  invisible(x)
}



#' @keywords internal
sigma.mo.bqr.svy <- function(x) {
  if (!inherits(x, "mo.bqr.svy")) stop("Not a 'mo.bqr.svy' object.", call. = FALSE)
  # prefer top-level x$sigma if present
  if (is.list(x$sigma) && length(x$sigma) == length(x$quantile)) return(x$sigma)
  # fallback: build from nested structure
  out <- vector("list", length(x$fit))
  names(out) <- paste0("tau=", formatC(x$quantile, digits = 3, format = "f"))
  for (qi in seq_along(x$fit)) {
    block <- x$fit[[qi]]
    if (!is.list(block$directions)) next
    vals <- vapply(block$directions, function(dk) {
      v <- tryCatch(as.numeric(dk$sigma)[1], error = function(e) NA_real_)
      ifelse(is.finite(v), v, NA_real_)
    }, numeric(1))
    names(vals) <- paste0("dir", seq_along(vals))
    out[[qi]] <- vals
  }
  out
}

#' @rdname print.bayesQRsurvey
#' @title Print a \code{mo.bqr.svy} model
#' @exportS3Method print mo.bqr.svy
print.mo.bqr.svy <- function(x, ...) {
  rule <- strrep("\u2500", 56)

  cat("\n")
  cat("  Multiple-Output Bayesian Quantile Regression\n")
  cat("  ", rule, "\n", sep = "")

  qs <- paste(formatC(x$quantile, digits = 3, format = "f"), collapse = ", ")
  cat("  Formula    : ", deparse(x$formula), "\n", sep = "")
  cat("  Quantiles  : ", qs, "\n", sep = "")
  cat("  Directions : ", x$n_dir, "\n", sep = "")
  if (!is.null(x$n_obs) && !is.null(x$response_dim))
    cat("  Sample     : ", x$n_obs, " obs, ", x$response_dim, " responses\n", sep = "")
  if (!is.null(x$estimate_sigma) && isFALSE(x$estimate_sigma))
    cat("  Scale      : sigma fixed at 1\n", sep = "")

  # Number of coefficients per direction
  is_new <- is.list(x$coefficients) &&
    all(vapply(x$coefficients, is.matrix, logical(1)))
  if (is_new) {
    n_coef <- nrow(x$coefficients[[1]])
    cat("  Coefficients: ", n_coef, " per direction\n", sep = "")
  }

  cat("  ", rule, "\n", sep = "")
  cat("  Use summary() for coefficients and convergence details.\n\n")

  invisible(x)
}

# Helper: print coefficient table transposed (directions as rows)
# so it fits in a standard console width
#' @noRd
.print_coef_table <- function(M, digits = 3, max_dir = 8, indent = 4) {
  n_dir   <- ncol(M)
  n_param <- nrow(M)
  pad     <- strrep(" ", indent)

  param_names <- rownames(M)
  if (is.null(param_names)) param_names <- paste0("V", seq_len(n_param))
  dir_names <- colnames(M)
  if (is.null(dir_names)) dir_names <- paste0("dir", seq_len(n_dir))

  # Transpose: directions as rows, parameters as columns
  Mt <- t(round(M, digits))
  vals <- format(Mt, digits = digits, nsmall = digits, trim = TRUE)

  # Column widths
  col_w <- pmax(nchar(param_names),
                apply(vals, 2, function(v) max(nchar(v))))
  dir_w <- max(nchar(dir_names), 3)

  # How many to show
  n_show <- min(n_dir, max_dir)

  # Header
  cat(pad, formatC("", width = dir_w), sep = "")
  for (j in seq_len(n_param))
    cat("  ", formatC(param_names[j], width = col_w[j], format = "s"), sep = "")
  cat("\n")

  # Rows
  for (i in seq_len(n_show)) {
    cat(pad, formatC(dir_names[i], width = dir_w, format = "s"), sep = "")
    for (j in seq_len(n_param))
      cat("  ", formatC(vals[i, j], width = col_w[j], format = "s"), sep = "")
    cat("\n")
  }
  if (n_dir > max_dir)
    cat(pad, "... ", n_dir - max_dir, " more directions ",
        "(use print(x, max_dir = ", n_dir, ") to show all)\n", sep = "")
}
