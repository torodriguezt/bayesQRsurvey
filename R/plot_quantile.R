#' Plot Method for Bayesian Weighted Quantile Regression (Survey)
#'
#' @description
#' Plot method for objects of class \code{bqr.svy} produced by \code{bqr.svy()}.
#' It can display fitted quantile curves, coefficient–quantile profiles,
#' MCMC trace plots, and posterior densities.
#'
#' @details
#' Supported plot types:
#' \itemize{
#'   \item \code{type = "fit"}: Fitted quantile curves versus a single
#'         numeric predictor. Optionally overlay observed points and credible
#'         bands. Other covariates can be held fixed via \code{at}.
#'   \item \code{type = "quantile"}: A single coefficient as a function
#'         of the quantile \eqn{\tau}. Optionally add a reference line at 0 and
#'         the corresponding OLS estimate.
#'   \item \code{type = "trace"}: MCMC trace for one selected
#'         coefficient at a chosen \eqn{\tau}.
#'   \item \code{type = "density"}: Posterior density for one selected
#'         coefficient at a chosen \eqn{\tau}.
#' }
#'
#' Notes:
#' \itemize{
#'   \item \code{tau} must be included in \code{x$quantile}. If \code{NULL}, all
#'         available quantiles in the object are used.
#'   \item For \code{type = "fit"}, \code{predictor} must be a numeric column in
#'         the original model. If \code{NULL}, the first numeric predictor
#'         (different from the response) is chosen automatically.
#'   \item For \code{type = "fit"}, \code{at} is a named list
#'         (\code{list(var = value, ...)}) used to fix other covariates while
#'         plotting versus \code{predictor}. Provide valid levels for factors.
#'   \item When \code{use_ggplot = TRUE}, a ggplot object is returned and the
#'         appearance is controlled by \code{theme_style} and
#'         \code{color_palette}. Otherwise, base graphics are used and the
#'         function returns \code{invisible(NULL)}.
#' }
#'
#' @param x Object of class \code{bqr.svy}.
#' @param y Ignored (S3 signature).
#' @param type One of \code{"fit"}, \code{"quantile"}, \code{"trace"},
#'   \code{"density"}.
#' @param predictor (fit) Name of a numeric predictor; if \code{NULL}, the first
#'   numeric predictor (excluding the response) is used.
#' @param tau Quantile(s) to plot; must appear in \code{x$quantile}. If
#'   \code{NULL}, all available are used.
#' @param which (quantile/trace/density) Coefficient name or index to display.
#' The default is the first coefficient associated with the first variable in the model.
#' @param add_points (fit) Logical; overlay observed data points.
#' @param combine (fit) Logical; if multiple \code{tau}: \code{TRUE} overlays
#'   curves in one panel; \code{FALSE} uses one panel per quantile.
#' @param show_ci (fit) Logical; draw credible bands.
#' @param ci_probs (fit) Length-2 numeric vector with lower/upper probabilities
#'   for credible bands.
#' @param at (fit) Named list of fixed values for non-\code{predictor}
#'   covariates (see Details).
#' @param grid_length (fit) Integer; number of points in the predictor grid.
#' @param points_alpha (fit) Point transparency in \code{[0,1]}.
#' @param point_size (fit) Point size.
#' @param line_size (fit/quantile) Line width for fitted/summary lines.
#' @param main Optional main title.
#' @param use_ggplot Logical; if \code{TRUE}, return a ggplot object.
#' @param theme_style (ggplot) One of \code{"minimal"}, \code{"classic"},
#'   \code{"bw"}, \code{"light"}.
#' @param color_palette (ggplot) One of \code{"viridis"}, \code{"plasma"},
#'   \code{"set2"}, \code{"dark2"}.
#' @param add_h0 (quantile) Logical; add a horizontal reference at \eqn{y = 0}.
#' @param add_ols (quantile) Logical; add the OLS estimate (dotted line) for the
#'   selected coefficient.
#' @param ols_fit (quantile) Optional precomputed \code{lm} object; if
#'   \code{NULL}, an \code{lm()} is fitted internally using \code{x$model} and
#'   \code{x$terms}.
#' @param ols_weights (quantile) Optional numeric vector of weights when fitting
#'   OLS internally (length must match \code{nrow(x$model)}).
#' @param ... Accepted for compatibility; ignored by internal plotting code.
#'
#' @return \code{invisible(NULL)} for base R graphics, or a ggplot object if
#'   \code{use_ggplot = TRUE}.
#'
#' @examples
#' \donttest{
#' data(mtcars)
#' fit <- bqr.svy(mpg ~ wt + hp + cyl, data = mtcars,
#'                quantile = c(0.5), method = "ald",
#'                niter = 20000, burnin = 10000, thin = 5)
#'
#' plot(fit, type = "fit", predictor = "wt", show_ci = TRUE)
#' plot(fit, type = "quantile", which = "wt", add_h0 = TRUE, add_ols = TRUE)
#' plot(fit, type = "trace", which = "wt", tau = 0.5)
#' plot(fit, type = "density", which = "wt", tau = 0.5)
#' }
#'
#' @aliases plot
#' @method plot bqr.svy
#' @rdname plot.bqr.svy
#' @name plot.bqr.svy
#' @export
plot.bqr.svy <- function(
    x, y = NULL,
    type = c("fit", "quantile", "trace", "density"),
    predictor = NULL,
    tau = NULL,
    which = NULL,
    add_points = TRUE,
    combine = TRUE,
    show_ci = FALSE,
    ci_probs = c(0.1, 0.9),
    at = NULL,
    grid_length = 200,
    points_alpha = 0.4,
    point_size = 1.5,
    line_size = 1.2,
    main = NULL,
    use_ggplot = TRUE,
    theme_style = c("minimal", "classic", "bw", "light"),
    color_palette = c("viridis", "plasma", "set2", "dark2"),
    add_h0 = FALSE,
    add_ols = FALSE,
    ols_fit = NULL,
    ols_weights = NULL,
    ...
) {
  type <- match.arg(type)
  theme_style <- match.arg(theme_style)
  color_palette <- match.arg(color_palette)

  is_multi <- inherits(x, "bwqr_fit_multi")
  taus_all <- as.numeric(x$quantile)
  if (is.null(tau)) tau <- taus_all
  tau <- sort(intersect(taus_all, unique(as.numeric(tau))))
  if (!length(tau)) stop("Requested 'tau' does not exist in the object.", call. = FALSE)

  mf <- x$model
  tt <- x$terms
  if (is.null(mf) || is.null(tt)) stop("Object does not contain 'model' and/or 'terms'.", call. = FALSE)

  X_colnames <- colnames(stats::model.matrix(tt, mf))
  resp <- as.character(stats::formula(tt))[2]

  # Helpers --------------------------
  .get_draws <- function(obj, tau_sel = NULL) {
    D <- if (inherits(obj, "bwqr_fit_multi")) {
      idx <- which.min(abs(obj$quantile - tau_sel))
      obj$draws[[idx]]
    } else obj$draws
    D <- as.matrix(D)
    keep <- intersect(colnames(D), X_colnames)
    if (!length(keep)) stop("The 'draws' do not contain the expected coefficient columns.", call. = FALSE)
    D[, keep, drop = FALSE]
  }
  .alpha <- function(col, a) grDevices::adjustcolor(col, alpha.f = max(min(a, 1), 0))
  .default_at <- function(df) {
    res <- list()
    for (nm in names(df)) {
      if (nm %in% c(resp)) next
      v <- df[[nm]]
      if (is.numeric(v))       res[[nm]] <- stats::median(v)
      else if (is.factor(v))   res[[nm]] <- levels(v)[1L]
      else if (is.logical(v))  res[[nm]] <- FALSE
      else                     res[[nm]] <- v[1L]
    }
    res
  }
  .make_newdata <- function(pred, at_vals, grid_len) {
    base_vals <- .default_at(mf)
    if (!is.null(at_vals)) base_vals[names(at_vals)] <- at_vals
    nd <- mf[rep(1, grid_len), , drop = FALSE]
    for (nm in names(base_vals)) {
      if (nm %in% names(nd) && nm != pred && nm != resp) nd[[nm]] <- base_vals[[nm]]
    }
    if (!is.numeric(mf[[pred]])) stop("The 'predictor' must be numeric.", call. = FALSE)
    xr <- range(mf[[pred]], na.rm = TRUE)
    nd[[pred]] <- seq(xr[1], xr[2], length.out = grid_len)
    nd
  }

  # Theme and color setup for ggplot
  .get_theme <- function(style) {
    switch(style,
           "minimal" = ggplot2::theme_minimal(),
           "classic" = ggplot2::theme_classic(),
           "bw" = ggplot2::theme_bw(),
           "light" = ggplot2::theme_light()
    )
  }

  .get_color_scale <- function(palette, n) {
    switch(palette,
           "viridis" = ggplot2::scale_color_viridis_d(option = "D"),
           "plasma" = ggplot2::scale_color_viridis_d(option = "C"),
           "set2" = ggplot2::scale_color_brewer(type = "qual", palette = "Set2"),
           "dark2" = ggplot2::scale_color_brewer(type = "qual", palette = "Dark2")
    )
  }

  .get_fill_scale <- function(palette) {
    switch(palette,
           "viridis" = ggplot2::scale_fill_viridis_d(option = "D", alpha = 0.3),
           "plasma" = ggplot2::scale_fill_viridis_d(option = "C", alpha = 0.3),
           "set2" = ggplot2::scale_fill_brewer(type = "qual", palette = "Set2", alpha = 0.3),
           "dark2" = ggplot2::scale_fill_brewer(type = "qual", palette = "Dark2", alpha = 0.3)
    )
  }

  # Colors (base R)
  cols <- if (exists("hcl.colors", where = asNamespace("grDevices"), inherits = FALSE)) {
    grDevices::hcl.colors(length(tau), "Dark 3")
  } else {
    grDevices::rainbow(length(tau))
  }

  if (type == "fit") {
    if (is.null(predictor)) {
      cand <- setdiff(names(mf), resp)
      predictor <- cand[which(vapply(mf[cand], is.numeric, logical(1)))[1]]
    }
    if (is.null(predictor) || !(predictor %in% names(mf)))
      stop("Could not determine 'predictor'. Pass it (predictor='...').", call. = FALSE)

    newdata <- .make_newdata(predictor, at, grid_length)
    Xg <- stats::model.matrix(tt, newdata)
    if (!all(colnames(Xg) %in% X_colnames)) {
      stop("The design matrix of 'newdata' does not match the fit.", call. = FALSE)
    }
    Xg <- Xg[, X_colnames, drop = FALSE]

    if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
      plot_data_list <- lapply(seq_along(tau), function(k) {
        ti <- tau[k]
        Dk <- .get_draws(x, tau_sel = ti)
        preds <- Xg %*% t(Dk)

        xg <- newdata[[predictor]]
        qmed <- apply(preds, 1, stats::median)

        df <- data.frame(
          x = xg,
          y = qmed,
          tau = sprintf("tau = %.3f", ti),
          tau_numeric = ti
        )

        if (show_ci) {
          qs <- t(apply(preds, 1, stats::quantile, probs = ci_probs))
          df$y_lower <- qs[, 1]
          df$y_upper <- qs[, 2]
        }
        df
      })

      plot_data <- do.call(rbind, plot_data_list)

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y))

      if (show_ci) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$y_lower, ymax = .data$y_upper,
                                                   fill = .data$tau), alpha = 0.3)
        p <- p + .get_fill_scale(color_palette)
      }

      if (add_points) {
        obs_data <- data.frame(x = mf[[predictor]], y = mf[[resp]])
        p <- p + ggplot2::geom_point(data = obs_data, ggplot2::aes(x = .data$x, y = .data$y),
                                     color = "gray30", alpha = points_alpha, size = point_size,
                                     inherit.aes = FALSE)
      }

      p <- p + ggplot2::geom_line(ggplot2::aes(color = .data$tau), linewidth = line_size)
      p <- p + .get_color_scale(color_palette, length(tau))

      p <- p + .get_theme(theme_style) +
        ggplot2::labs(
          x = predictor,
          y = "Predicted Quantile",
          color = "Quantile",
          fill = "Quantile",
          title = if (is.null(main)) {
            if (length(tau) == 1L) sprintf("Quantile Regression Fit (tau = %.3f)", tau[1])
            else "Quantile Regression Fits"
          } else main
        )

      if (length(tau) > 1 && !combine) {
        p <- p + ggplot2::facet_wrap(~ .data$tau, scales = "free_y")
      }

      p <- p + ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        legend.position = "bottom",
        legend.box = "horizontal",
        panel.grid.minor = ggplot2::element_blank()
      )
      return(p)

    } else {
      return(.plot_fit_base(x, predictor, tau, mf, tt, X_colnames, resp, newdata, Xg,
                            add_points, combine, show_ci, ci_probs, grid_length,
                            points_alpha, point_size, line_size, main, cols))
    }
  }

  if (type == "quantile") {
    if (length(tau) < 2L) {
      stop("For 'type=\"quantile\"' you must have at least two quantiles in the object or pass 'tau' with length > 1.", call. = FALSE)
    }
    # Coefficient selection
    D_example <- .get_draws(x, tau_sel = tau[1])
    if (is.null(which)) which <- colnames(D_example)[1]
    idx <- if (is.character(which)) match(which[1], colnames(D_example)) else which[1]
    if (is.na(idx) || idx < 1 || idx > ncol(D_example)) stop("'which' out of range.", call. = FALSE)
    coef_name <- colnames(D_example)[idx]

    # Summary by quantile
    qsum <- lapply(tau, function(ti) {
      Dk <- .get_draws(x, tau_sel = ti)[, idx]
      c(tau = ti,
        med = stats::median(Dk),
        lo  = stats::quantile(Dk, probs = ci_probs[1]),
        hi  = stats::quantile(Dk, probs = ci_probs[2]))
    })
    qsum <- as.data.frame(do.call(rbind, qsum))

    # Optional OLS
    ols_coef <- NA_real_
    if (isTRUE(add_ols)) {
      if (is.null(ols_fit)) {
        fm <- stats::formula(tt)
        if (is.null(ols_weights)) ols_fit <- stats::lm(fm, data = mf)
        else                      ols_fit <- stats::lm(fm, data = mf, weights = ols_weights)
      }
      cn <- names(stats::coef(ols_fit))
      j  <- if (is.character(which)) match(coef_name, cn) else idx
      if (!is.na(j) && j >= 1 && j <= length(cn)) ols_coef <- stats::coef(ols_fit)[j]
    }

    # ggplot2 plot
    if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(qsum, ggplot2::aes(x = tau, y = med)) +
        ggplot2::geom_line(size = line_size) +
        ggplot2::geom_point(size = 2)

      if (isTRUE(show_ci)) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                                      alpha = 0.25, inherit.aes = TRUE)
      }
      if (isTRUE(add_h0))  p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "solid")
      if (isTRUE(add_ols) && is.finite(ols_coef)) {
        p <- p + ggplot2::geom_hline(yintercept = ols_coef, linetype = "dotted")
      }

      p <- p + .get_theme(theme_style) +
        ggplot2::labs(
          x = "Quantile",
          y = coef_name,
          title = if (is.null(main)) sprintf("Coefficient across quantiles: %s", coef_name) else main
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "none"
        )
      return(p)

    } else {
      # Base R
      ylim <- range(qsum$lo, qsum$hi, qsum$med, if (add_h0) 0 else NA_real_, na.rm = TRUE)
      graphics::plot(qsum$tau, qsum$med, type = "n",
                     xlab = "Quantile", ylab = coef_name,
                     main = if (is.null(main)) sprintf("Coefficient across quantiles: %s", coef_name) else main,
                     ylim = ylim)
      graphics::grid()
      if (isTRUE(show_ci)) {
        xx <- c(qsum$tau, rev(qsum$tau))
        yy <- c(qsum$lo,  rev(qsum$hi))
        graphics::polygon(xx, yy, col = grDevices::adjustcolor("gray50", 0.3), border = NA)
      }
      graphics::lines(qsum$tau, qsum$med, lwd = line_size)
      graphics::points(qsum$tau, qsum$med, pch = 19)
      if (isTRUE(add_h0))  graphics::abline(h = 0, lty = 1)
      if (isTRUE(add_ols) && is.finite(ols_coef)) graphics::abline(h = ols_coef, lty = 3)
      invisible(NULL)
    }
  }


  if (type == "trace") {
    D <- .get_draws(x, tau_sel = tau[1])
    if (is.null(which)) which <- colnames(D)[1]
    idx <- if (is.character(which)) match(which[1], colnames(D)) else which[1]
    if (is.na(idx) || idx < 1 || idx > ncol(D)) stop("'which' out of range.", call. = FALSE)

    if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
      v <- D[, idx]
      nm <- colnames(D)[idx]
      trace_data <- data.frame(
        iteration = seq_along(v),
        value = v
      )

      p <- ggplot2::ggplot(trace_data, ggplot2::aes(x = .data$iteration, y = .data$value))
      p <- p + ggplot2::geom_line(color = "steelblue", size = 0.7, alpha = 0.8)
      p <- p + ggplot2::geom_hline(yintercept = stats::median(v), color = "red",
                                   linetype = "dashed", size = 1)

      p <- p + .get_theme(theme_style)
      p <- p + ggplot2::labs(
        x = "Iteration",
        y = nm,
        title = if (is.null(main)) sprintf("MCMC Trace: %s (tau = %.3f)", nm, tau[1]) else main
      )

      p <- p + ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10),
        panel.grid.minor = ggplot2::element_blank()
      )
      return(p)

    } else {
      return(.plot_trace_base(D, idx, tau, main))
    }
  }

  if (type == "density") {
    D <- .get_draws(x, tau_sel = tau[1])
    if (is.null(which)) which <- colnames(D)[1]
    idx <- if (is.character(which)) match(which[1], colnames(D)) else which[1]
    if (is.na(idx) || idx < 1 || idx > ncol(D)) stop("'which' out of range.", call. = FALSE)

    if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
      v <- D[, idx]
      nm <- colnames(D)[idx]

      p <- ggplot2::ggplot(data.frame(x = v), ggplot2::aes(x = .data$x))
      p <- p + ggplot2::geom_density(fill = "lightblue", color = "darkblue",
                                     alpha = 0.7, size = 1)
      p <- p + ggplot2::geom_vline(xintercept = stats::median(v), color = "red",
                                   linetype = "dashed", size = 1)

      p <- p + .get_theme(theme_style)
      p <- p + ggplot2::labs(
        x = nm,
        y = "Density",
        title = if (is.null(main)) sprintf("Posterior Density: %s (tau = %.3f)", nm, tau[1]) else main
      )

      p <- p + ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10),
        panel.grid.minor = ggplot2::element_blank()
      )
      return(p)

    } else {
      return(.plot_density_base(D, idx, tau, main))
    }
  }

  invisible(NULL)
}


# --------------------------------------------------------------------
# Helper functions para base R
# --------------------------------------------------------------------

.plot_fit_base <- function(x, predictor, tau, mf, tt, X_colnames, resp, newdata, Xg,
                           add_points, combine, show_ci, ci_probs, grid_length,
                           points_alpha, point_size, line_size, main, cols) {

  .get_draws <- function(obj, tau_sel = NULL) {
    D <- if (inherits(obj, "bwqr_fit_multi")) {
      idx <- which.min(abs(obj$quantile - tau_sel))
      obj$draws[[idx]]
    } else obj$draws
    D <- as.matrix(D)
    keep <- intersect(colnames(D), X_colnames)
    D[, keep, drop = FALSE]
  }

  .alpha <- function(col, a) grDevices::adjustcolor(col, alpha.f = max(min(a, 1), 0))
  .best_mfrow <- function(n) { r <- floor(sqrt(n)); c <- ceiling(n / r); c(r, c) }

  preds_list <- lapply(tau, function(ti) {
    Dk <- .get_draws(x, tau_sel = ti)
    Xg %*% t(Dk)
  })
  y_rng <- range(do.call(cbind, preds_list), na.rm = TRUE)
  xg <- newdata[[predictor]]

  if (length(tau) == 1L || isTRUE(combine)) {
    graphics::plot(xg, apply(preds_list[[1]], 1, stats::median), type = "n",
                   xlab = predictor, ylab = "Predicted quantile",
                   main = if (is.null(main)) {
                     if (length(tau) == 1L) sprintf("Quantile fit vs %s (tau=%.3f)", predictor, tau[1])
                     else sprintf("Quantile fit vs %s (taus: %s)", predictor, paste(format(tau, digits = 3), collapse = ", "))
                   } else main,
                   ylim = y_rng)
    graphics::grid()

    if (isTRUE(show_ci)) {
      q_low <- ci_probs[1]; q_high <- ci_probs[2]
      for (k in seq_along(tau)) {
        preds_k <- preds_list[[k]]
        y_low <- apply(preds_k, 1, stats::quantile, probs = q_low)
        y_high <- apply(preds_k, 1, stats::quantile, probs = q_high)
        polygon_col <- .alpha(cols[k], 0.3)
        graphics::polygon(c(xg, rev(xg)), c(y_low, rev(y_high)),
                          col = polygon_col, border = NA)
      }
    }

    for (k in seq_along(tau)) {
      preds_k <- preds_list[[k]]
      y_med <- apply(preds_k, 1, stats::median)
      graphics::lines(xg, y_med, col = cols[k], lwd = line_size)
    }

    if (isTRUE(add_points)) {
      pred_vals <- mf[[predictor]]
      resp_vals <- mf[[resp]]
      pts_col <- .alpha("gray30", points_alpha)
      graphics::points(pred_vals, resp_vals, col = pts_col, pch = 19, cex = point_size)
    }

    if (length(tau) > 1L) {
      labs <- sprintf("tau = %.3f", tau)
      graphics::legend("topright", legend = labs, col = cols, lwd = line_size, bty = "n")
    }

  } else {
    mfr <- .best_mfrow(length(tau))
    op <- graphics::par(mfrow = mfr)
    on.exit(graphics::par(op), add = TRUE)

    for (k in seq_along(tau)) {
      ti <- tau[k]
      preds_k <- preds_list[[k]]
      y_med <- apply(preds_k, 1, stats::median)

      graphics::plot(xg, y_med, type = "l", col = cols[k], lwd = line_size,
                     xlab = predictor, ylab = "Predicted quantile",
                     main = sprintf("tau = %.3f", ti))
      graphics::grid()

      if (isTRUE(show_ci)) {
        q_low <- ci_probs[1]; q_high <- ci_probs[2]
        y_low <- apply(preds_k, 1, stats::quantile, probs = q_low)
        y_high <- apply(preds_k, 1, stats::quantile, probs = q_high)
        polygon_col <- .alpha(cols[k], 0.3)
        graphics::polygon(c(xg, rev(xg)), c(y_low, rev(y_high)),
                          col = polygon_col, border = NA)
      }

      if (isTRUE(add_points)) {
        pred_vals <- mf[[predictor]]
        resp_vals <- mf[[resp]]
        pts_col <- .alpha("gray30", points_alpha)
        graphics::points(pred_vals, resp_vals, col = pts_col, pch = 19, cex = point_size)
      }
    }
  }
  invisible(NULL)
}

.plot_trace_base <- function(D, idx, tau, main) {
  v <- D[, idx]
  nm <- colnames(D)[idx]
  graphics::plot(v, type = "l", col = "steelblue",
                 main = if (is.null(main)) sprintf("MCMC trace: %s (tau=%.3f)", nm, tau[1]) else main,
                 xlab = "Iteration", ylab = nm)
  graphics::abline(h = stats::median(v), col = "red", lty = 2)
  graphics::grid(nx = NA, ny = NULL)
  invisible(NULL)
}

.plot_density_base <- function(D, idx, tau, main) {
  v <- D[, idx]
  nm <- colnames(D)[idx]
  d <- stats::density(v)
  graphics::plot(d, main = if (is.null(main)) sprintf("Posterior density: %s (tau=%.3f)", nm, tau[1]) else main,
                 xlab = nm)
  graphics::abline(v = stats::median(v), col = "red", lty = 2)
  graphics::grid(nx = NA, ny = NULL)
  invisible(NULL)
}


#' @rdname plot.bqr.svy
#' @export
plot.bwqr_fit <- function(x, ...) {
  plot.bqr.svy(x, ...)
}

#' @rdname plot.bqr.svy
#' @export
plot.bwqr_fit_multi <- function(x, ...) {
  plot.bqr.svy(x, ...)
}
