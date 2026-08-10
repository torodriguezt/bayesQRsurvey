## Replication script for
##   "bayesQRsurvey: Bayesian Quantile Regression for Complex Survey Data"
##
## Reproduces every figure, output and table in the manuscript.
## Figures are written to ./Figures. Runtime: about two minutes.

for (lib in .libPaths()) {
  if (dir.exists(file.path(lib, "bayesQRsurvey")))
    remove.packages("bayesQRsurvey", lib = lib)
}

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

remotes::install_github("torodriguezt/bayesQRsurvey",
                        upgrade = "never", force = TRUE)

library("bayesQRsurvey")
library("ggplot2")
library("patchwork")


## Cosmetic settings, not shown in the manuscript. They affect only the
## appearance of the saved figures, as do the ggsave() calls below.

OUT <- "Figures"
dir.create(OUT, showWarnings = FALSE)

theme_single <- theme_classic(base_size = 22) +
  theme(axis.text              = element_text(colour = "black"),
        legend.title           = element_blank(),
        legend.text            = element_text(size = 18),
        legend.position        = "inside",
        legend.position.inside = c(0.83, 0.17),
        legend.background      = element_blank())

theme_region <- theme_classic(base_size = 14) +
  theme(axis.text              = element_text(colour = "black"),
        legend.position        = "inside",
        legend.position.inside = c(0.14, 0.86),
        legend.background      = element_blank(),
        legend.key.height      = unit(1.4, "lines"))

update_geom_defaults("point", list(size = 3))

axis_small <- theme(axis.title = element_text(size = 14),
                    axis.text  = element_text(size = 12))
axis_wide  <- theme(axis.title = element_text(size = 20),
                    axis.text  = element_text(size = 17))

diag_style <- function(p) p & theme_classic(base_size = 13) &
  theme(axis.text = element_text(colour = "black"))


## Section 4: colour palettes used by both examples.

scale_sex <- scale_colour_manual(values = c(Boys = "black",
                                            Girls = "grey70"))

scale_tau <- scale_colour_manual(
  values = c("0.1" = "grey70", "0.5" = "grey40", "0.9" = "black"),
  breaks = c("0.9", "0.5", "0.1"),
  labels = c("0.9" = expression(tau == 0.9),
             "0.5" = expression(tau == 0.5),
             "0.1" = expression(tau == 0.1)))


## Section 4.1: Example 1, single-output.

theme_set(theme_single)

## The Anthro data, with age in years and sex labelled.

data("Anthro", package = "bayesQRsurvey")
str(Anthro)

Anthro$age <- Anthro$age / 12
Anthro$sex <- factor(Anthro$sex, levels = c("1", "0"),
                     labels = c("Boys", "Girls"))

## Figure 1: weight and height against age, by sex.

p1 <- ggplot(Anthro, aes(x = age, y = wgt, colour = sex)) +
  geom_point() + scale_sex +
  labs(x = "Age (years)", y = "Weight (kg)")

p2 <- ggplot(Anthro, aes(x = age, y = hgt, colour = sex)) +
  geom_point() + scale_sex +
  labs(x = "Age (years)", y = "Height (cm)")

p1 + p2

ggsave(file.path(OUT, "figure2.pdf"),
       (p1 + axis_wide) + (p2 + axis_wide), width = 12.5, height = 5)

## Fit with the asymmetric Laplace method.

set.seed(50)
fit_ald <- bqr.svy(wgt ~ age + I(age^2) + sex, weights = dweight,
                   data = Anthro, quantile = c(0.1, 0.5, 0.9),
                   niter = 20000, burnin = 10000, thin = 1)

fit_ald
print(summary(fit_ald), tau = 0.5)
fit_ald$diagnosis[["tau=0.500"]]

## Figure 2: trace plots (top) and posterior densities (bottom) at tau = 0.5.

p4 <- plot(fit_ald, type = "trace", tau = 0.5,
           color_palette = "grey", theme_style = "none")

p5 <- plot(fit_ald, type = "density", tau = 0.5,
           color_palette = "grey", theme_style = "none")

p4
p5

ggsave(file.path(OUT, "figure4_trace.pdf"),   diag_style(p4),
       width = 9, height = 6)
ggsave(file.path(OUT, "figure5_density.pdf"), diag_style(p5),
       width = 9, height = 6)

## Figure 3: fitted quantile curves for weight against age.

p3 <- plot(fit_ald, type = "fit", which = "age", add_points = FALSE,
           color_palette = "none", theme_style = "none") +
  scale_tau + labs(x = "Age (years)", y = "Weight (kg)")

p3

ggsave(file.path(OUT, "figure3.pdf"),
       p3 + axis_small + theme(legend.key.height = unit(1.5, "lines")),
       width = 6.25, height = 5)

## Figure 4: coefficients across a grid of quantiles, with the OLS estimate.

quantvec <- seq(0.1, 0.9, by = 0.05)

set.seed(50)
fit_ald_grid <- bqr.svy(wgt ~ age + I(age^2) + sex, weights = dweight,
                        data = Anthro, quantile = quantvec,
                        niter = 50000, burnin = 25000, thin = 1)

p10 <- plot(fit_ald_grid, type = "quantile",
            which = c("(Intercept)", "age", "I(age^2)", "sexGirls"),
            add_ols = TRUE, color_palette = "grey", theme_style = "none") +
  labs(x = "quantile")

p10

ggsave(file.path(OUT, "figure10_quantile_ald.pdf"), diag_style(p10),
       width = 8, height = 6)

## The same model under the two remaining methods.

set.seed(50)
fit_score <- bqr.svy(wgt ~ age + I(age^2) + sex,
                     weights = dweight, data = Anthro, method = "score",
                     quantile = c(0.1, 0.5, 0.9),
                     niter = 50000, burnin = 10000, thin = 1)

set.seed(50)
fit_ap <- bqr.svy(wgt ~ age + I(age^2) + sex, weights = dweight,
                  data = Anthro, method = "approximate",
                  quantile = c(0.1, 0.5, 0.9),
                  niter = 350000, burnin = 50000, thin = 50)

## Fit under an informative prior.

myprior <- prior(beta_x_mean = rep(0, 4), beta_x_cov = 25)

set.seed(50)
fit_ald_prior <- bqr.svy(wgt ~ age + I(age^2) + sex, weights = dweight,
                         data = Anthro, quantile = c(0.1, 0.5, 0.9),
                         niter = 20000, burnin = 10000, thin = 1,
                         prior = myprior)


## Section 4.2: Example 2, multiple-output.

theme_set(theme_region)

## Bivariate fit over a grid of 20 directions.

set.seed(50)
fit_mo <- mo.bqr.svy(cbind(wgt, hgt) ~ age + I(age^2) + sex,
                     weights = dweight, data = Anthro,
                     quantile = c(0.05, 0.10, 0.15),
                     n_dir = 20, max_iter = 2000)

fit_mo
print(summary(fit_mo), coefficients = FALSE)

## Figure 5: nested quantile regions for a two-year-old boy.
## plotQuantileRegion() returns a list, so the plot is taken from $plot.

reg <- plotQuantileRegion(fit_mo, response = c("wgt", "hgt"),
                          datafile = Anthro, xValue = c(1, 2, 4, 0),
                          ngridpoints = 450, paintedArea = FALSE,
                          color_palette = "grey", theme_style = "none")

reg$plot

ggsave(file.path(OUT, "plotQuantileRegion1.pdf"), reg$plot,
       width = 6.5, height = 5.5)

## The same fit with the directions supplied through U and gamma_U.

n_dir   <- 20
angles  <- (0:(n_dir - 1)) * 2 * pi / n_dir
U       <- rbind(cos(angles), sin(angles))
gamma_U <- lapply(seq_len(n_dir), function(k)
  matrix(c(-sin(angles[k]), cos(angles[k])), ncol = 1))

set.seed(50)
fit_mo_manual <- mo.bqr.svy(cbind(wgt, hgt) ~ age + I(age^2) + sex,
                            weights = dweight, data = Anthro,
                            quantile = c(0.05, 0.10, 0.15),
                            U = U, gamma_U = gamma_U, max_iter = 2000)


## Appendix, Table 1: the three methods side by side.

summary(fit_ald)
summary(fit_score)
summary(fit_ap)


## Computational environment.

sessionInfo()
