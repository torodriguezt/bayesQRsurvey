#' Plot method for Bayesian quantile regression fits
#'
#' @description
#' Plot method for objects of class \code{bqr.svy} produced by \code{bqr.svy()}.
#' It can display fitted quantile curves, coefficient-quantile profiles,
#' MCMC trace plots, and posterior densities.
#'
#' @details
#' The \code{type} argument selects one of four displays.
#' \itemize{
#'   \item \code{"fit"} draws the fitted quantile curves against a numeric
#'         predictor, optionally with the observed data and credible bands.
#'   \item \code{"quantile"} draws one or more coefficients as functions of
#'         the quantile \eqn{\tau}, with a credible band.
#'   \item \code{"trace"} draws the MCMC traces at one quantile.
#'   \item \code{"density"} draws the posterior densities at one quantile.
#' }
#'
#' @param x an object of class \code{"bqr.svy"}.
#' @param type one of \code{"fit"}, \code{"quantile"}, \code{"trace"} or
#'   \code{"density"} (see Details).
#' @param tau quantile(s) to plot, among those fitted. If \code{NULL}, all
#'   fitted quantiles are used.
#' @param which what to display. For \code{type = "fit"}, the name of the
#'   numeric predictor on the horizontal axis, by default the first one in the
#'   model. For the other types, coefficient names or indices, by default the
#'   first coefficient for \code{"quantile"} and all of them for
#'   \code{"trace"} and \code{"density"}.
#' @param add_points logical; overlay the observed data on the fitted curves.
#' @param combine logical; with several quantiles, draw the fitted curves in one
#'   panel (\code{TRUE}) or in one panel per quantile (\code{FALSE}).
#' @param show_ci logical; add credible bands to the fitted curves.
#' @param ci_probs lower and upper probabilities of the credible intervals
#'   drawn in the plots. Default \code{c(0.025, 0.975)}.
#' @param at named list of values for the covariates not on the horizontal axis
#'   of the fit plot, for example \code{list(sex = "Girls")}. Covariates not
#'   listed are held at their median, or at their first level for factors.
#' @param grid_length number of predictor values at which the fitted curves are
#'   evaluated.
#' @param main optional plot title.
#' @param theme_style \code{"minimal"} (default) or \code{"none"}. With
#'   \code{"none"} no theme is applied, so the active \code{theme_set()} is used.
#' @param color_palette one of \code{"viridis"}, \code{"plasma"},
#'   \code{"set2"}, \code{"dark2"}, \code{"grey"} for greyscale printing, or
#'   \code{"none"} to leave the colour scale to the user.
#' @param add_h0 logical; add a horizontal line at zero to the quantile plot.
#' @param add_ols logical; add the least-squares estimate to the quantile plot.
#' @param ols_weights optional weights for the least-squares fit drawn by
#'   \code{add_ols}, one per observation. By default the fit is unweighted.
#' @param ... not used; any argument passed here is an error.
#'
#' @return A ggplot object.
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
#' @method plot bqr.svy
#' @rdname plot.bqr.svy
#' @name plot.bqr.svy
#' @export
plot.bqr.svy <- function(
    x,
    type = c("fit", "quantile", "trace", "density"),
    tau = NULL,
    which = NULL,
    add_points = TRUE,
    combine = TRUE,
    show_ci = FALSE,
    ci_probs = c(0.025, 0.975),
    at = NULL,
    grid_length = 200,
    main = NULL,
    theme_style = c("minimal", "none"),
    color_palette = c("viridis", "plasma", "set2", "dark2", "grey", "none"),
    add_h0 = FALSE,
    add_ols = FALSE,
    ols_weights = NULL,
    ...
) {
  .stop_on_dots(match.call(expand.dots = FALSE)$...)
  type <- match.arg(type)
  theme_style <- match.arg(theme_style)
  color_palette <- match.arg(color_palette)

  taus_all <- as.numeric(x$quantile)
  if (is.null(tau)) tau <- taus_all
  tau <- sort(intersect(taus_all, unique(as.numeric(tau))))
  if (!length(tau)) stop("Requested 'tau' does not exist in the object.", call. = FALSE)

  mf <- x$model
  tt <- x$terms
  if (is.null(mf) || is.null(tt)) stop("Object does not contain 'model' and/or 'terms'.", call. = FALSE)

  X_colnames <- colnames(.model_matrix(x))
  resp <- as.character(stats::formula(tt))[2]

  # Helpers --------------------------
  .get_draws <- function(obj, tau_sel = NULL) {
    D <- .draws_for_tau(obj, tau_sel, include_sigma = FALSE)
    keep <- intersect(colnames(D), X_colnames)
    if (!length(keep)) stop("The 'draws' do not contain the expected coefficient columns.", call. = FALSE)
    D[, keep, drop = FALSE]
  }
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

  # "none" applies no theme, so the plot honours the active theme_set() (useful
  # for reproducing a paper's own styling). "minimal" is a light preset that
  # inherits the base theme's colours and fonts rather than forcing greys.
  .get_theme <- function(style) {
    if (identical(style, "none")) return(NULL)
    ggplot2::theme_minimal() + ggplot2::theme(
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
    ref_col      <- "black"    # reference lines (posterior mean, OLS): black + dashed/dotted
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

  if (type == "fit") {
    predictor <- which
    if (is.null(predictor)) {
      cand <- setdiff(names(mf), resp)
      predictor <- cand[vapply(mf[cand], is.numeric, logical(1))][1]
    }
    if (is.null(predictor) || !(predictor %in% names(mf)))
      stop("Could not determine predictor. Pass it via 'which'.", call. = FALSE)

    newdata <- .make_newdata(predictor, at, grid_length)
    Xg <- stats::model.matrix(stats::delete.response(tt), newdata,
                              contrasts.arg = x$contrasts)
    if (!all(colnames(Xg) %in% X_colnames)) {
      stop("The design matrix of 'newdata' does not match the fit.", call. = FALSE)
    }
    Xg <- Xg[, X_colnames, drop = FALSE]

    plot_data_list <- lapply(seq_along(tau), function(k) {
      ti <- tau[k]
      Dk <- .get_draws(x, tau_sel = ti)
      preds <- Xg %*% t(Dk)

      xg <- newdata[[predictor]]
      # the curve is the conditional quantile at the posterior mean, x' beta_bar,
      # the same quantity fitted() returns at the observed covariates
      qfit <- drop(Xg %*% colMeans(Dk))

      df <- data.frame(
        x = xg,
        y = qfit,
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
                                   color = "grey50", alpha = 0.4, size = 2.5,
                                   shape = 16, inherit.aes = FALSE)
    }

    p <- p + ggplot2::geom_line(ggplot2::aes(color = .data$tau), linewidth = 1.2)
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

    # Optional least-squares reference. Fitted directly on the design matrix of
    # the model, so the weights are taken as given rather than looked up through
    # lm()'s non-standard evaluation.
    ols_coef <- NULL
    if (isTRUE(add_ols)) {
      X_ols <- .model_matrix(x)
      y_ols <- stats::model.response(mf)
      if (is.null(ols_weights)) {
        ols_coef <- stats::lm.fit(X_ols, y_ols)$coefficients
      } else {
        if (!is.numeric(ols_weights) || length(ols_weights) != nrow(X_ols))
          stop("'ols_weights' must be a numeric vector with one weight per observation (",
               nrow(X_ols), ").", call. = FALSE)
        ols_coef <- stats::lm.wfit(X_ols, y_ols, ols_weights)$coefficients
      }
    }

    # Summary by quantile for each coefficient: the posterior mean, which is
    # what coef() and summary() report, with a credible band from the draws
    qsum_list <- lapply(seq_along(which_idx), function(w) {
      idx <- which_idx[w]
      coef_name <- which_names[w]
      do.call(rbind, lapply(tau, function(ti) {
        Dk <- .get_draws(x, tau_sel = ti)[, idx]
        data.frame(
          tau = ti,
          est = mean(Dk),
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
    if (!is.null(ols_coef)) {
      ols_vals <- unname(ols_coef[which_names])
      ols_df <- data.frame(coef = factor(which_names, levels = which_names),
                           ols = ols_vals, stringsAsFactors = FALSE)
      ols_df <- ols_df[is.finite(ols_df$ols), , drop = FALSE]
    }

    multi_coef <- length(which_idx) > 1L

    p <- ggplot2::ggplot(qsum, ggplot2::aes(x = .data$tau, y = .data$est))

    # Always show the credible band for quantile plots
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                                  fill = accent_fill, alpha = 0.3,
                                  inherit.aes = TRUE)

    if (isTRUE(add_h0)) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                   color = "grey50", linewidth = 0.5)
    }
    if (!is.null(ols_df) && nrow(ols_df) > 0L) {
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
  }

  # type is "trace" or "density": both show the draws at a single quantile
  D <- .get_draws(x, tau_sel = tau[1])
  # which = NULL selects every coefficient; a name/index vector selects a subset
  if (is.null(which)) which <- colnames(D)
  idx <- if (is.character(which)) match(which, colnames(D)) else which
  if (anyNA(idx) || any(idx < 1) || any(idx > ncol(D)))
    stop("'which' out of range.", call. = FALSE)
  nms <- colnames(D)[idx]
  # reference line at the posterior mean, the estimate coef() reports
  est_data <- data.frame(coef = factor(nms, levels = nms),
                         est  = colMeans(D[, idx, drop = FALSE]))

  if (type == "trace") {
    trace_data <- do.call(rbind, lapply(seq_along(idx), function(j)
      data.frame(iteration = seq_len(nrow(D)), value = D[, idx[j]],
                 coef = nms[j], stringsAsFactors = FALSE)))
    trace_data$coef <- factor(trace_data$coef, levels = nms)

    p <- ggplot2::ggplot(trace_data,
                         ggplot2::aes(x = .data$iteration, y = .data$value))
    p <- p + ggplot2::geom_line(color = accent_col, linewidth = 0.3, alpha = 0.7)
    p <- p + ggplot2::geom_hline(data = est_data,
                                 ggplot2::aes(yintercept = .data$est),
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
  }

  dens_data <- do.call(rbind, lapply(seq_along(idx), function(j)
    data.frame(x = D[, idx[j]], coef = nms[j], stringsAsFactors = FALSE)))
  dens_data$coef <- factor(dens_data$coef, levels = nms)
  ci_data <- do.call(rbind, lapply(seq_along(idx), function(j)
    data.frame(coef = nms[j],
               xint = as.numeric(stats::quantile(D[, idx[j]], probs = ci_probs)),
               stringsAsFactors = FALSE)))
  ci_data$coef <- factor(ci_data$coef, levels = nms)

  p <- ggplot2::ggplot(dens_data, ggplot2::aes(x = .data$x))
  p <- p + ggplot2::geom_density(fill = accent_fill, color = accent_col,
                                 alpha = 0.6, linewidth = 0.7)
  p <- p + ggplot2::geom_vline(data = est_data,
                               ggplot2::aes(xintercept = .data$est),
                               color = ref_col, linetype = "dashed",
                               linewidth = 0.8)
  p <- p + ggplot2::geom_vline(data = ci_data,
                               ggplot2::aes(xintercept = .data$xint),
                               color = "grey50", linetype = "dotted",
                               linewidth = 0.5)
  if (length(idx) > 1L)
    p <- p + ggplot2::facet_wrap(~ .data$coef, scales = "free")
  p <- p + .get_theme(theme_style)
  p + ggplot2::labs(
    x = if (length(idx) == 1L) nms[1] else "Value",
    y = "Density",
    title = main
  )
}


# Removed or renamed plot arguments would otherwise fall into `...` and be
# ignored without notice, changing the figure. Fail as R does for an unused
# argument, naming the replacement when there is one.
.stop_on_dots <- function(dots, renamed = character()) {
  if (!length(dots)) return(invisible(NULL))
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))
  hint <- ifelse(nms %in% names(renamed),
                 paste0(" (now '", renamed[nms], "')"), "")
  shown <- ifelse(nzchar(nms), paste0("'", nms, "'"), "an unnamed argument")
  stop("Unused argument(s): ", paste0(shown, hint, collapse = ", "), ".",
       call. = FALSE)
}
