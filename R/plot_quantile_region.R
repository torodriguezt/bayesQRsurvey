# ---- internal helpers ------------------------------------------------

#' Check which grid points fall inside the quantile region
#' @keywords internal
#' @noRd
.checkPoints_R <- function(seqY1, seqY2, directions, orthBases,
                           betas, xValue) {
  n_dirs <- ncol(directions)
  p      <- length(xValue)

  a_list <- matrix(0, nrow = 2, ncol = n_dirs)
  b_vals <- numeric(n_dirs)

  for (d in seq_len(n_dirs)) {
    u    <- directions[, d]
    v    <- orthBases[, d]
    beta <- betas[, d]
    a_list[, d] <- u - beta[p + 1] * v
    b_vals[d]   <- sum(xValue * beta[seq_len(p)])
  }

  grid  <- expand.grid(y1 = seqY1, y2 = seqY2)
  y_mat <- as.matrix(grid)
  lhs   <- y_mat %*% a_list

  inside <- rep(TRUE, nrow(grid))
  for (d in seq_len(n_dirs)) inside <- inside & (lhs[, d] >= b_vals[d])

  if (sum(inside) > 0L) as.matrix(grid[inside, ]) else matrix(nrow = 0L, ncol = 2L)
}

#' Extract beta coefficient matrix for a given quantile index
#' @keywords internal
#' @noRd
.extract_betas_mobqr <- function(model, tau_idx) {
  tau_name <- names(model$fit)[tau_idx]
  model$coefficients[[tau_name]]
}

# ---- main function ---------------------------------------------------

#' Plot Bivariate Quantile Regions for Multiple-Output Models
#'
#' @description
#' Draws bivariate quantile regions from a fitted \code{mo.bqr.svy} object.
#' The function projects data onto a grid, determines which points lie inside
#' each quantile region, and visualises the result using \pkg{ggplot2}.
#'
#' @details
#' Two display modes are available:
#' \itemize{
#'    \item \code{paintedArea = TRUE} (default): Filled ribbons with contour
#'          outlines, using a sequential \code{"Blues"} palette.  An in-plot
#'          text legend shows the quantile level and approximate coverage.
#'    \item \code{paintedArea = FALSE}: Contour-only lines coloured by
#'          quantile, with a standard ggplot2 legend at the bottom.
#' }
#'
#' The quantile regions are built from the directional approach described in
#' Nascimento & \enc{Gonçalves}{Goncalves} (2025).  For each quantile
#' \eqn{\tau}, the function evaluates every grid point against the
#' half-space constraints defined by the fitted directions and orthogonal
#' bases.
#'
#' @param model An object of class \code{"mo.bqr.svy"} produced by
#'    \code{\link{mo.bqr.svy}}.
#' @param response Character vector of length 2 giving the names of the
#'    response columns in \code{datafile}.
#' @param datafile A data frame containing at least the two columns named
#'    in \code{response}.
#' @param ngridpoints Integer; number of grid points per axis
#'    (default = 200).
#' @param xValue Numeric vector of covariate values at which to evaluate
#'    the regression.  For an intercept-only model use \code{xValue = 1}.
#' @param paintedArea Logical; if \code{TRUE} (default) the regions are
#'    drawn as filled ribbons (see Details).
#' @param range_y An optional \eqn{2 \times 2} matrix whose rows give
#'    the \code{[min, max]} range for the first and second response,
#'    respectively.  If \code{NULL}, the ranges are computed from the data
#'    with a 15 percent margin.
#' @param main Optional character string for the plot title.
#' @param theme_style Character; one of \code{"minimal"}, \code{"classic"},
#'    \code{"bw"}, \code{"light"} (default = \code{"bw"}).
#' @param color_palette Character; one of \code{"blues"}, \code{"set1"},
#'    \code{"viridis"}, \code{"dark2"} (default = \code{"blues"}).
#' @param point_alpha Numeric in \code{[0, 1]}; transparency for the
#'    observed-data point cloud (default = 0.3).
#' @param point_size Numeric; size of the data points (default = 1.2).
#' @param line_size Numeric; line width for contour outlines
#'    (default = 0.8).
#' @param verbose Logical; if \code{TRUE}, print per-quantile progress
#'    messages (default = \code{FALSE}).
#'
#' @return Invisibly, a list with components:
#'    \item{plot}{A \code{ggplot} object.}
#'    \item{data}{A data frame with columns \code{y1}, \code{min},
#'      \code{max} and \code{tau}, one block per quantile level.}
#'
#' @examples
#' \donttest{
#' # Load the Anthro dataset (children anthropometric data, already cleaned)
#' data(Anthro, package = "bayesQRsurvey")
#'
#' # --- Directions ---
#' n_dir  <- 12
#' angles <- (0:(n_dir - 1)) * 2 * pi / n_dir
#' U_mat  <- rbind(cos(angles), sin(angles))
#' gamma_list <- lapply(seq_len(n_dir), function(k)
#'    matrix(c(-sin(angles[k]), cos(angles[k])), ncol = 1))
#'
#' # --- Fit ---
#' fit <- mo.bqr.svy(
#'    cbind(wgt, hgt) ~ 1,
#'    data     = Anthro,
#'    quantile = c(0.05, 0.10, 0.25),
#'    U        = U_mat,
#'    gamma_U  = gamma_list,
#'    max_iter = 2000,
#'    verbose  = FALSE
#' )
#'
#' # --- Plot quantile regions (filled ribbons) ---
#' plotQuantileRegion(fit, response = c("wgt", "hgt"), datafile = Anthro)
#'
#' # Contour-only style
#' plotQuantileRegion(fit, response = c("wgt", "hgt"), datafile = Anthro,
#'                    paintedArea = FALSE)
#' }
#'
#' @seealso \code{\link{mo.bqr.svy}}, \code{\link{prior}}
#' @importFrom ggplot2 ggplot aes geom_point geom_ribbon geom_path
#'    scale_color_brewer scale_linetype_discrete theme theme_bw theme_minimal
#'    theme_classic theme_light labs annotate element_text element_blank
#'    scale_fill_brewer scale_color_viridis_d scale_fill_viridis_d
#' @export
plotQuantileRegion <- function(model,
                               response,
                               datafile,
                               ngridpoints   = 200,
                               xValue        = 1,
                               paintedArea   = TRUE,
                               range_y       = NULL,
                               main          = NULL,
                               theme_style   = c("bw", "minimal", "classic", "light"),
                               color_palette = c("blues", "set1", "viridis", "dark2"),
                               point_alpha   = 0.3,
                               point_size    = 1.2,
                               line_size     = 0.8,
                               verbose       = FALSE) {

  # --- Validate inputs ------------------------------------------------
  if (!inherits(model, "mo.bqr.svy"))
    stop("'model' must be of class 'mo.bqr.svy'.", call. = FALSE)
  if (!is.character(response) || length(response) != 2L)
    stop("'response' must be a character vector of length 2.", call. = FALSE)
  if (!is.data.frame(datafile))
    stop("'datafile' must be a data frame.", call. = FALSE)
  miss <- setdiff(response, names(datafile))
  if (length(miss))
    stop(sprintf("Column(s) %s not found in 'datafile'.",
                 paste0("'", miss, "'", collapse = ", ")), call. = FALSE)

  theme_style   <- match.arg(theme_style)
  color_palette <- match.arg(color_palette)

  # --- Extract model components ---------------------------------------
  directions <- model$U
  K          <- ncol(directions)
  d          <- nrow(directions)

  orthBases <- matrix(0, nrow = d, ncol = K)
  for (k in seq_len(K)) {
    gk <- model$Gamma_list[[k]]
    if (is.matrix(gk) && ncol(gk) >= 1L) {
      orthBases[, k] <- gk[, 1]
    } else if (is.numeric(gk) && length(gk) == d) {
      orthBases[, k] <- gk
    }
  }

  taus  <- model$quantile
  ntaus <- length(taus)

  betaDifDirections <- lapply(seq_len(ntaus), function(qi) {
    .extract_betas_mobqr(model, qi)
  })

  # --- Observed response matrix (complete cases) ----------------------
  Y <- datafile[, response, drop = FALSE]
  Y <- Y[stats::complete.cases(Y), , drop = FALSE]

  if (is.null(range_y)) {
    y1range <- range(Y[[1]], na.rm = TRUE)
    y2range <- range(Y[[2]], na.rm = TRUE)
    y1range <- y1range + c(-1, 1) * diff(y1range) * 0.15
    y2range <- y2range + c(-1, 1) * diff(y2range) * 0.15
  } else {
    y1range <- range_y[1, ]
    y2range <- range_y[2, ]
  }

  seqY1 <- seq(y1range[1], y1range[2], length.out = ngridpoints)
  seqY2 <- seq(y2range[1], y2range[2], length.out = ngridpoints)

  # --- Compute regions per tau ----------------------------------------
  regions <- list()
  for (idx in seq_len(ntaus)) {
    betas_mat <- betaDifDirections[[idx]]
    if (any(is.na(betas_mat))) {
      if (verbose)
        message(sprintf("  tau = %.2f: skipped (NA coefficients)", taus[idx]))
      next
    }

    cp <- .checkPoints_R(seqY1, seqY2, directions, orthBases,
                         betas_mat, xValue)

    if (nrow(cp) > 0L) {
      y1_in  <- sort(unique(cp[, 1]))
      y2_max <- vapply(y1_in, function(v) max(cp[cp[, 1] == v, 2]),
                       numeric(1))
      y2_min <- vapply(y1_in, function(v) min(cp[cp[, 1] == v, 2]),
                       numeric(1))
      regions[[as.character(taus[idx])]] <- data.frame(
        y1 = y1_in, min = y2_min, max = y2_max, tau = taus[idx]
      )
      if (verbose)
        message(sprintf("  tau = %.2f: %d grid slices", taus[idx],
                        length(y1_in)))
    } else {
      if (verbose)
        message(sprintf("  tau = %.2f: empty region", taus[idx]))
    }
  }

  all_regions <- do.call(rbind, regions)
  if (is.null(all_regions) || nrow(all_regions) == 0L) {
    warning("No quantile regions could be computed.", call. = FALSE)
    return(invisible(NULL))
  }

  # --- ggplot2 theme helper -------------------------------------------
  base_theme <- switch(theme_style,
                       "minimal" = ggplot2::theme_minimal(),
                       "classic" = ggplot2::theme_classic(),
                       "bw"      = ggplot2::theme_bw(),
                       "light"   = ggplot2::theme_light()
  )

  # --- Observed data points -------------------------------------------
  df_points <- data.frame(y1 = Y[[1]], y2 = Y[[2]])

  # --- Build plot -----------------------------------------------------
  if (isTRUE(paintedArea)) {
    # --- Painted-area style -------------------------------------------
    taus_order <- sort(unique(all_regions$tau))
    n_taus     <- length(taus_order)

    # Sequential blue palette via grDevices (no RColorBrewer dependency)
    # Skip the 2 lightest shades so all colours remain readable
    n_extra <- 2L
    all_cols <- grDevices::hcl.colors(n_taus + n_extra, "Blues 3")
    palette  <- all_cols[(n_extra + 1L):(n_taus + n_extra)]

    # Darker palette for legend text so labels are readable on white bg
    legend_cols <- grDevices::hcl.colors(n_taus + n_extra, "Blues 3",
                                         alpha = NULL)
    # Darken by shifting luminance down
    legend_palette <- vapply(palette, function(col) {
      rgb_val <- grDevices::col2rgb(col)
      grDevices::rgb(
        pmax(rgb_val[1, ] - 60L, 0L),
        pmax(rgb_val[2, ] - 60L, 0L),
        pmax(rgb_val[3, ] - 60L, 0L),
        maxColorValue = 255
      )
    }, character(1), USE.NAMES = FALSE)

    g <- ggplot2::ggplot(df_points, ggplot2::aes(x = .data$y1, y = .data$y2)) +
      ggplot2::geom_point(alpha = point_alpha, color = "gray40",
                          size = point_size) +
      base_theme

    for (i in seq_along(taus_order)) {
      reg <- all_regions[all_regions$tau == taus_order[i], , drop = FALSE]

      g <- g +
        ggplot2::geom_ribbon(
          data        = reg,
          ggplot2::aes(x = .data$y1, ymin = .data$min, ymax = .data$max),
          alpha       = 0.3,
          fill        = palette[i],
          inherit.aes = FALSE
        ) +
        ggplot2::geom_path(
          data = data.frame(
            x = c(reg$y1, rev(reg$y1), reg$y1[1]),
            y = c(reg$min, rev(reg$max), reg$min[1])
          ),
          ggplot2::aes(x = .data$x, y = .data$y),
          color       = palette[i],
          linewidth   = line_size,
          inherit.aes = FALSE
        )
    }

    legend_labels <- paste0("tau = ", taus_order,
                            " (~", round((1 - taus_order) * 100), "%)")

    g <- g +
      ggplot2::annotate(
        "text",
        x        = y1range[1] + diff(y1range) * 0.02,
        y        = y2range[2] - diff(y2range) * 0.04 * seq_len(n_taus),
        label    = legend_labels,
        color    = legend_palette,
        hjust    = 0,
        size     = 3.5,
        fontface = "bold"
      ) +
      ggplot2::labs(
        x     = response[1],
        y     = response[2],
        title = if (is.null(main)) "Bivariate Quantile Regions" else main
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 12),
        axis.text  = ggplot2::element_text(size = 10)
      )

  } else {
    # --- Contour-only style -------------------------------------------
    contour_data <- do.call(rbind, lapply(unique(all_regions$tau), function(tv) {
      reg <- all_regions[all_regions$tau == tv, , drop = FALSE]
      data.frame(
        x    = c(reg$y1, rev(reg$y1), reg$y1[1]),
        y    = c(reg$min, rev(reg$max), reg$min[1]),
        taus = tv
      )
    }))

    .contour_color_scale <- switch(color_palette,
                                   "blues"   = ggplot2::scale_color_brewer(palette = "Blues",
                                                                           name = expression(tau)),
                                   "set1"    = ggplot2::scale_color_brewer(palette = "Set1",
                                                                           name = expression(tau)),
                                   "viridis" = ggplot2::scale_color_viridis_d(option = "D",
                                                                              name = expression(tau)),
                                   "dark2"   = ggplot2::scale_color_brewer(palette = "Dark2",
                                                                           name = expression(tau))
    )

    g <- ggplot2::ggplot() +
      base_theme +
      ggplot2::geom_point(
        data  = df_points,
        ggplot2::aes(x = .data$y1, y = .data$y2),
        alpha = point_alpha,
        color = "gray50",
        size  = point_size
      ) +
      ggplot2::geom_path(
        data = contour_data,
        ggplot2::aes(x = .data$x, y = .data$y,
                     color    = factor(.data$taus),
                     linetype = factor(.data$taus)),
        linewidth = line_size
      ) +
      .contour_color_scale +
      ggplot2::scale_linetype_discrete(name = expression(tau)) +
      ggplot2::labs(
        x     = response[1],
        y     = response[2],
        title = if (is.null(main)) "Bivariate Quantile Regions" else main
      ) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title      = ggplot2::element_text(size = 12),
        axis.text       = ggplot2::element_text(size = 10),
        legend.position = "bottom",
        legend.box      = "vertical"
      )
  }

  print(g)
  invisible(list(plot = g, data = all_regions))
}
