#' Plot Method for Bayesian Weighted Quantile Regression
#'
#' @description
#' Plot method for objects of class \code{bqr.svy} produced by \code{bqr.svy()}.
#' It can display fitted quantile curves, coefficient-quantile profiles,
#' MCMC trace plots, and posterior densities.
#'
#' @details
#' Supported plot types:
#' \itemize{
#'   \item \code{type = "fit"}: Fitted quantile curves versus a single
#'         numeric predictor (selected via \code{which}). Optionally overlay
#'         observed points and credible bands. Other covariates can be held
#'         fixed via \code{at}.
#'   \item \code{type = "quantile"}: A single coefficient as a function
#'         of the quantile \eqn{\tau}. Optionally add a reference line at 0 and
#'         the corresponding OLS estimate.
#'   \item \code{type = "trace"}: MCMC trace(s) at a chosen \eqn{\tau}. One or
#'         several coefficients may be selected via \code{which}; when more
#'         than one is shown the panels are arranged with \code{facet_wrap}.
#'   \item \code{type = "density"}: Posterior density(ies) at a chosen
#'         \eqn{\tau}, faceted over the selected coefficients as for
#'         \code{type = "trace"}.
#' }
#'
#' Notes:
#' \itemize{
#'   \item \code{tau} must be included in \code{x$quantile}. If \code{NULL}, all
#'         available quantiles in the object are used.
#'   \item For \code{type = "fit"}, \code{which} must name a numeric column in
#'         the original model. If \code{NULL}, the first numeric predictor
#'         (different from the response) is chosen automatically.
#'   \item For \code{type = "fit"}, \code{at} is a named list
#'         (\code{list(var = value, ...)}) used to fix other covariates while
#'         plotting versus the selected predictor. Provide valid levels for
#'         factors.
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
#' @param tau Quantile(s) to plot; must appear in \code{x$quantile}. If
#'   \code{NULL}, all available are used.
#' @param which Variable name(s) or coefficient index(es) to display.
#'   For \code{type = "fit"}, the name of a numeric predictor to plot on the
#'   x-axis (if \code{NULL}, the first numeric predictor is used).
#'   For \code{type = "quantile"}, a character vector of coefficient names
#'   or integer vector of indices; when more than one is given the plot uses
#'   \code{facet_wrap} to show all coefficients in a single figure.
#'   For \code{type = "trace"} and \code{type = "density"}, a character vector
#'   of coefficient names or integer vector of indices; when \code{NULL} (the
#'   default) all coefficients are shown, and when more than one is selected the
#'   panels are arranged with \code{facet_wrap}.
#' @param add_points (fit) Logical; overlay observed data points.
#' @param combine (fit) Logical; if multiple \code{tau}: \code{TRUE} overlays
#'   curves in one panel; \code{FALSE} uses one panel per quantile.
#' @param show_ci (fit) Logical; draw credible bands.
#' @param ci_probs Length-2 numeric vector with the lower/upper probabilities of
#'   the credible interval shown by the fit ribbon, the density plot bounds, and
#'   the quantile-process band. Default \code{c(0.025, 0.975)} (95\%).
#' @param at (fit) Named list of fixed values for non-plotted
#'   covariates (see Details).
#' @param grid_length (fit) Integer; number of points in the predictor grid.
#' @param points_alpha (fit) Point transparency in \code{[0,1]}.
#' @param point_size (fit) Point size.
#' @param line_size (fit/quantile) Line width for fitted/summary lines.
#' @param main Optional main title.
#' @param use_ggplot Logical; if \code{TRUE}, return a ggplot object.
#' @param theme_style (ggplot) One of \code{"minimal"}, \code{"classic"},
#'   \code{"bw"}, \code{"light"}, or \code{"none"}. With \code{"none"} the plot
#'   applies no theme and therefore honours the active \code{theme_set()}.
#' @param color_palette (ggplot) One of \code{"viridis"}, \code{"plasma"},
#'   \code{"set2"}, \code{"dark2"}, \code{"grey"} for a print-friendly greyscale
#'   rendering (also switches the single-colour accents used by the trace,
#'   density, and quantile plots to greys), or \code{"none"} to add no colour
#'   scale so a user-supplied \code{scale_colour_*()} controls the mapping.
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
#'                quantile = c(0.5, 0.75), method = "ald",
#'                niter = 20000, burnin = 10000, thin = 5)
#'
#' plot(fit, type = "fit", which = "wt", show_ci = TRUE)
#' plot(fit, type = "quantile", which = "wt", add_h0 = TRUE, add_ols = TRUE)
#' plot(fit, type = "quantile", which = c("(Intercept)", "wt", "hp", "cyl"),
#'      add_h0 = TRUE, add_ols = TRUE)
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
    tau = NULL,
    which = NULL,
    add_points = TRUE,
    combine = TRUE,
    show_ci = FALSE,
    ci_probs = c(0.025, 0.975),
    at = NULL,
    grid_length = 200,
    points_alpha = 0.4,
    point_size = 2.5,
    line_size = 1.2,
    main = NULL,
    use_ggplot = TRUE,
    theme_style = c("minimal", "classic", "bw", "light", "none"),
    color_palette = c("viridis", "plasma", "set2", "dark2", "grey", "none"),
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
  .make_newdata <- function(pred, at_vals, grid_len) {
    # Build the prediction grid from the base variables of the model formula,
    # not from the (already transformed) model-frame columns. The model frame
    # keeps its "terms" attribute even after subsetting, so reusing its columns
    # makes model.matrix() treat the grid as a ready-made model frame and freeze
    # transformed terms such as I(age^2) or log(x). Rebuilding from the base
    # variables lets model.matrix() re-evaluate those terms over the grid, so
    # curved fits are drawn correctly.
    base_vars <- intersect(all.vars(stats::delete.response(tt)), names(mf))
    if (!(pred %in% base_vars))
      stop("Predictor '", pred, "' is not a base variable of the model.",
           call. = FALSE)
    if (!is.numeric(mf[[pred]]))
      stop("The 'predictor' must be numeric.", call. = FALSE)

    xr <- range(mf[[pred]], na.rm = TRUE)
    nd <- list()
    for (nm in base_vars) {
      v <- mf[[nm]]
      if (nm == pred) {
        nd[[nm]] <- seq(xr[1], xr[2], length.out = grid_len)
      } else if (!is.null(at_vals) && !is.null(at_vals[[nm]])) {
        nd[[nm]] <- if (is.factor(v)) factor(at_vals[[nm]], levels = levels(v))
                    else at_vals[[nm]]
      } else if (is.factor(v)) {
        nd[[nm]] <- factor(levels(v)[1L], levels = levels(v))
      } else if (is.logical(v)) {
        nd[[nm]] <- FALSE
      } else if (is.numeric(v)) {
        nd[[nm]] <- stats::median(v, na.rm = TRUE)
      } else {
        nd[[nm]] <- v[1L]
      }
    }
    as.data.frame(nd, stringsAsFactors = FALSE)
  }

  # Theme setup for ggplot. "none" applies no theme, so the plot honours the
  # active theme_set() (useful for reproducing a paper's own styling). The
  # other styles use a light preset that inherits the base theme's colours and
  # fonts rather than forcing greys, so the output stays clean and composable.
  .get_theme <- function(style) {
    if (identical(style, "none")) return(NULL)
    base <- switch(style,
           "minimal" = ggplot2::theme_minimal(),
           "classic" = ggplot2::theme_classic(),
           "bw" = ggplot2::theme_bw(),
           "light" = ggplot2::theme_light()
    )
    base + ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(12, 12, 8, 8)
    )
  }

  # Single-color accent for 1-tau plots; multi-tau uses the palette.
  # The "grey" palette swaps these accents for print-friendly greys (the
  # reference lines stay distinguishable via their dashed/dotted linetypes).
  if (identical(color_palette, "grey")) {
    accent_col   <- "grey20"   # near-black for lines/points and density borders
    accent_fill  <- "grey70"   # light grey for ribbons/density fills
    ref_col      <- "black"    # reference lines (median, OLS): black + dashed/dotted
  } else {
    accent_col   <- "#2171B5"   # strong blue
    accent_fill  <- "#6BAED6"   # lighter blue
    ref_col      <- "#D6604D"   # muted red for reference lines
  }

  .get_color_scale <- function(palette, n) {
    if (identical(palette, "none")) return(NULL)  # user supplies their own scale
    if (n == 1L) {
      return(ggplot2::scale_color_manual(values = accent_col))
    }
    switch(palette,
           "viridis" = ggplot2::scale_color_viridis_d(option = "D"),
           "plasma" = ggplot2::scale_color_viridis_d(option = "C"),
           "set2" = ggplot2::scale_color_brewer(type = "qual", palette = "Set2"),
           "dark2" = ggplot2::scale_color_brewer(type = "qual", palette = "Dark2"),
           "grey" = ggplot2::scale_color_grey(start = 0.7, end = 0.1)
    )
  }

  .get_fill_scale <- function(palette, n = 2L) {
    if (identical(palette, "none")) return(NULL)  # user supplies their own scale
    if (n == 1L) {
      return(ggplot2::scale_fill_manual(values = accent_fill))
    }
    switch(palette,
           "viridis" = ggplot2::scale_fill_viridis_d(option = "D", alpha = 0.3),
           "plasma" = ggplot2::scale_fill_viridis_d(option = "C", alpha = 0.3),
           "set2" = ggplot2::scale_fill_brewer(type = "qual", palette = "Set2", alpha = 0.3),
           "dark2" = ggplot2::scale_fill_brewer(type = "qual", palette = "Dark2", alpha = 0.3),
           "grey" = ggplot2::scale_fill_grey(start = 0.8, end = 0.3)
    )
  }

  # Colors (base R)
  cols <- if (identical(color_palette, "grey")) {
    grDevices::grey.colors(length(tau), start = 0.7, end = 0.1)
  } else if (exists("hcl.colors", where = asNamespace("grDevices"), inherits = FALSE)) {
    grDevices::hcl.colors(length(tau), "Dark 3")
  } else {
    grDevices::rainbow(length(tau))
  }

  if (type == "fit") {
    predictor <- which
    if (is.null(predictor)) {
      cand <- setdiff(names(mf), resp)
      predictor <- cand[vapply(mf[cand], is.numeric, logical(1))][1]
    }
    if (is.null(predictor) || !(predictor %in% names(mf)))
      stop("Could not determine predictor. Pass it via 'which'.", call. = FALSE)

    newdata <- .make_newdata(predictor, at, grid_length)
    Xg <- stats::model.matrix(stats::delete.response(tt), newdata)
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
          tau = format(ti, trim = TRUE),
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
                                                   fill = .data$tau), alpha = 0.25)
        p <- p + .get_fill_scale(color_palette, length(tau))
      }

      if (add_points) {
        obs_data <- data.frame(x = mf[[predictor]], y = mf[[resp]])
        p <- p + ggplot2::geom_point(data = obs_data, ggplot2::aes(x = .data$x, y = .data$y),
                                     color = "grey50", alpha = points_alpha, size = point_size,
                                     shape = 16, inherit.aes = FALSE)
      }

      p <- p + ggplot2::geom_line(ggplot2::aes(color = .data$tau), linewidth = line_size)
      p <- p + .get_color_scale(color_palette, length(tau))

      p <- p + .get_theme(theme_style) +
        ggplot2::labs(
          x = predictor,
          y = resp,
          color = expression(tau),
          title = main
        )

      if (show_ci) {
        p <- p + ggplot2::labs(fill = expression(tau))
      }

      if (length(tau) > 1 && !combine) {
        p <- p + ggplot2::facet_wrap(~ .data$tau, scales = "free_y")
      }

      # Legend placement is only imposed when the plot supplies its own theme;
      # with theme_style = "none" the active theme_set() controls the legend.
      if (theme_style != "none") {
        if (length(tau) > 1L) {
          p <- p + ggplot2::theme(
            legend.title = ggplot2::element_text(size = 11, face = "bold"),
            legend.text = ggplot2::element_text(size = 10),
            legend.position = "bottom",
            legend.box = "horizontal"
          )
        } else {
          p <- p + ggplot2::theme(legend.position = "none")
        }
      }
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
    # Coefficient selection - allow multiple coefficients
    D_example <- .get_draws(x, tau_sel = tau[1])
    if (is.null(which)) which <- colnames(D_example)[1]

    # Resolve indices for all requested coefficients
    which_names <- character(length(which))
    which_idx   <- integer(length(which))
    for (w in seq_along(which)) {
      idx_w <- if (is.character(which[w])) match(which[w], colnames(D_example)) else which[w]
      if (is.na(idx_w) || idx_w < 1 || idx_w > ncol(D_example))
        stop(sprintf("'which' value '%s' out of range.", which[w]), call. = FALSE)
      which_names[w] <- colnames(D_example)[idx_w]
      which_idx[w]   <- idx_w
    }

    # Optional OLS (compute once)
    ols_fit_obj <- NULL
    if (isTRUE(add_ols)) {
      if (is.null(ols_fit)) {
        fm <- stats::formula(tt)
        if (is.null(ols_weights)) ols_fit_obj <- stats::lm(fm, data = mf)
        else                      ols_fit_obj <- stats::lm(fm, data = mf, weights = ols_weights)
      } else {
        ols_fit_obj <- ols_fit
      }
    }

    # Summary by quantile for each coefficient
    qsum_list <- lapply(seq_along(which_idx), function(w) {
      idx <- which_idx[w]
      coef_name <- which_names[w]
      do.call(rbind, lapply(tau, function(ti) {
        Dk <- .get_draws(x, tau_sel = ti)[, idx]
        data.frame(
          tau = ti,
          med = stats::median(Dk),
          lo  = unname(stats::quantile(Dk, probs = ci_probs[1])),
          hi  = unname(stats::quantile(Dk, probs = ci_probs[2])),
          coef = coef_name,
          stringsAsFactors = FALSE
        )
      }))
    })
    qsum <- do.call(rbind, qsum_list)
    qsum$coef <- factor(qsum$coef, levels = which_names)

    # OLS reference values per coefficient
    ols_df <- NULL
    if (isTRUE(add_ols) && !is.null(ols_fit_obj)) {
      cn <- names(stats::coef(ols_fit_obj))
      ols_vals <- vapply(which_names, function(nm) {
        j <- match(nm, cn)
        if (!is.na(j)) stats::coef(ols_fit_obj)[j] else NA_real_
      }, numeric(1))
      ols_df <- data.frame(coef = factor(which_names, levels = which_names),
                           ols = ols_vals, stringsAsFactors = FALSE)
      ols_df <- ols_df[is.finite(ols_df$ols), , drop = FALSE]
    }

    # Single coefficient case (backwards compatible)
    multi_coef <- length(which_idx) > 1L

    # ggplot2 plot
    if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(qsum, ggplot2::aes(x = .data$tau, y = .data$med))

      # Always show the credible band for quantile plots
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                                    fill = accent_fill, alpha = 0.3,
                                    inherit.aes = TRUE)

      if (isTRUE(add_h0)) {
        p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                     color = "grey50", linewidth = 0.5)
      }
      if (isTRUE(add_ols) && !is.null(ols_df) && nrow(ols_df) > 0L) {
        p <- p + ggplot2::geom_hline(data = ols_df,
                                     ggplot2::aes(yintercept = .data$ols),
                                     linetype = "dotted",
                                     color = ref_col, linewidth = 0.8)
      }

      p <- p + ggplot2::geom_line(color = accent_col, linewidth = 0.5,
                                  linetype = "dashed")
      p <- p + ggplot2::geom_point(color = accent_col, size = 1.6, shape = 1)

      if (multi_coef) {
        p <- p + ggplot2::facet_wrap(~ .data$coef, scales = "free_y")
      }

      p <- p + .get_theme(theme_style) +
        ggplot2::labs(
          x = expression(tau),
          y = if (!multi_coef) which_names[1] else NULL,
          title = main
        ) +
        ggplot2::theme(legend.position = "none")

      # Bold facet strips are only imposed when the plot supplies its own theme;
      # with theme_style = "none" the active theme_set() controls the strips.
      if (multi_coef && theme_style != "none") {
        p <- p + ggplot2::theme(
          strip.text = ggplot2::element_text(size = 11, face = "bold")
        )
      }
      return(p)

    } else {
      # Base R (single coefficient only)
      idx <- which_idx[1]
      coef_name <- which_names[1]
      ols_coef <- if (!is.null(ols_df) && nrow(ols_df) > 0L) ols_df$ols[1] else NA_real_
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
    # which = NULL selects every coefficient; a name/index vector selects a subset
    if (is.null(which)) which <- colnames(D)
    idx <- if (is.character(which)) match(which, colnames(D)) else which
    if (anyNA(idx) || any(idx < 1) || any(idx > ncol(D)))
      stop("'which' out of range.", call. = FALSE)
    nms <- colnames(D)[idx]

    if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
      trace_data <- do.call(rbind, lapply(seq_along(idx), function(j)
        data.frame(iteration = seq_len(nrow(D)), value = D[, idx[j]],
                   coef = nms[j], stringsAsFactors = FALSE)))
      trace_data$coef <- factor(trace_data$coef, levels = nms)
      med_data <- data.frame(coef = factor(nms, levels = nms),
                             med  = apply(D[, idx, drop = FALSE], 2, stats::median))

      p <- ggplot2::ggplot(trace_data,
                           ggplot2::aes(x = .data$iteration, y = .data$value))
      p <- p + ggplot2::geom_line(color = accent_col, linewidth = 0.3, alpha = 0.7)
      p <- p + ggplot2::geom_hline(data = med_data,
                                   ggplot2::aes(yintercept = .data$med),
                                   color = ref_col, linetype = "dashed",
                                   linewidth = 0.8)
      if (length(idx) > 1L)
        p <- p + ggplot2::facet_wrap(~ .data$coef, scales = "free_y")
      p <- p + .get_theme(theme_style)
      p <- p + ggplot2::labs(
        x = "Iteration",
        y = if (length(idx) == 1L) nms[1] else "Value",
        title = main
      )
      return(p)

    } else {
      return(.plot_trace_base(D, idx, tau, main))
    }
  }

  if (type == "density") {
    D <- .get_draws(x, tau_sel = tau[1])
    # which = NULL selects every coefficient; a name/index vector selects a subset
    if (is.null(which)) which <- colnames(D)
    idx <- if (is.character(which)) match(which, colnames(D)) else which
    if (anyNA(idx) || any(idx < 1) || any(idx > ncol(D)))
      stop("'which' out of range.", call. = FALSE)
    nms <- colnames(D)[idx]

    if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
      dens_data <- do.call(rbind, lapply(seq_along(idx), function(j)
        data.frame(x = D[, idx[j]], coef = nms[j], stringsAsFactors = FALSE)))
      dens_data$coef <- factor(dens_data$coef, levels = nms)
      med_data <- data.frame(coef = factor(nms, levels = nms),
                             med  = apply(D[, idx, drop = FALSE], 2, stats::median))
      ci_data <- do.call(rbind, lapply(seq_along(idx), function(j)
        data.frame(coef = nms[j],
                   xint = as.numeric(stats::quantile(D[, idx[j]], probs = ci_probs)),
                   stringsAsFactors = FALSE)))
      ci_data$coef <- factor(ci_data$coef, levels = nms)

      p <- ggplot2::ggplot(dens_data, ggplot2::aes(x = .data$x))
      p <- p + ggplot2::geom_density(fill = accent_fill, color = accent_col,
                                     alpha = 0.6, linewidth = 0.7)
      p <- p + ggplot2::geom_vline(data = med_data,
                                   ggplot2::aes(xintercept = .data$med),
                                   color = ref_col, linetype = "dashed",
                                   linewidth = 0.8)
      p <- p + ggplot2::geom_vline(data = ci_data,
                                   ggplot2::aes(xintercept = .data$xint),
                                   color = "grey50", linetype = "dotted",
                                   linewidth = 0.5)
      if (length(idx) > 1L)
        p <- p + ggplot2::facet_wrap(~ .data$coef, scales = "free")
      p <- p + .get_theme(theme_style)
      p <- p + ggplot2::labs(
        x = if (length(idx) == 1L) nms[1] else "Value",
        y = "Density",
        title = main
      )
      return(p)

    } else {
      return(.plot_density_base(D, idx, tau, main))
    }
  }

  invisible(NULL)
}


# --------------------------------------------------------------------
# Helper functions for base R
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
                   xlab = predictor, ylab = resp,
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
                     xlab = predictor, ylab = resp,
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
  if (length(idx) > 1L) {
    op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
    r <- floor(sqrt(length(idx)))
    graphics::par(mfrow = c(r, ceiling(length(idx) / r)))
  }
  for (i in idx) {
    v <- D[, i]; nm <- colnames(D)[i]
    graphics::plot(v, type = "l", col = "steelblue",
                   main = if (is.null(main)) sprintf("MCMC trace: %s (tau=%.3f)", nm, tau[1]) else main,
                   xlab = "Iteration", ylab = nm)
    graphics::abline(h = stats::median(v), col = "red", lty = 2)
    graphics::grid(nx = NA, ny = NULL)
  }
  invisible(NULL)
}

.plot_density_base <- function(D, idx, tau, main) {
  if (length(idx) > 1L) {
    op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
    r <- floor(sqrt(length(idx)))
    graphics::par(mfrow = c(r, ceiling(length(idx) / r)))
  }
  for (i in idx) {
    v <- D[, i]; nm <- colnames(D)[i]
    d <- stats::density(v)
    graphics::plot(d, main = if (is.null(main)) sprintf("Posterior density: %s (tau=%.3f)", nm, tau[1]) else main,
                   xlab = nm)
    graphics::abline(v = stats::median(v), col = "red", lty = 2)
    graphics::grid(nx = NA, ny = NULL)
  }
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
