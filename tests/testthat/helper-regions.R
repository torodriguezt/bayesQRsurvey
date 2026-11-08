# The quantile regions drawn by plot.mo.bqr.svy: the data of its ribbon layer
# (painted_area = TRUE), one row per grid slice with columns y1, min, max, tau.
region_data <- function(p) {
  is_ribbon <- vapply(p$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))
  p$layers[[which(is_ribbon)]]$data
}
