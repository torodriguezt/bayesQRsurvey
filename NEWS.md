# bayesQRsurvey 0.2.1

## Bug fixes

* Fixed `plot()` crashing with "contrasts can be applied only to factors with 2 or more levels" when the model included factor predictors.
* Fixed DOI formatting in documentation (`\url{doi:...}` changed to `\doi{...}`), which caused help pages to fail to render.
* Fixed `Anthro.rda` compression to prevent `R_decompress1` warnings during installation.

## Data

* `Anthro$sex` is now stored as a factor (levels: `0` = girl, `1` = boy), so users no longer need to call `as.factor()` before fitting models.

## Plot improvements

* Refined color palette across all plot types: cohesive blue accent (`#2171B5`) for lines, lighter blue (`#6BAED6`) for fills, muted red (`#D6604D`) for reference lines.
* Single-quantile fit plots now use a solid blue line instead of a single viridis color, and hide the redundant legend.
* Quantile-coefficient plots (`type = "quantile"`) now always show the credible band and use the Greek tau symbol on the x-axis.
* Posterior density plots now display credible interval bounds as dotted vertical lines.
* Trace plots use thinner lines for better readability with dense MCMC chains.
* Legend labels for multiple quantiles now show only the numeric value (e.g., `0.100`) with the legend title as the tau symbol, instead of repeating "tau = " on every entry.
* Titles are now left-aligned with a subtitle for parameter details (variable name, tau value) in a lighter style.
* Removed grey panel background for a cleaner look.

# bayesQRsurvey 0.2.0

* Added `plotQuantileRegion` for multivariate quantile regression.
* Redesigned `print()` and `summary()` methods for better readability.
* Improved aesthetics and parameter handling in plotting functions.
* Internal cleanup of compiled code and documentation updates.
