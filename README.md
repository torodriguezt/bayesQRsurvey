# bayesQRsurvey

[![CRAN](https://img.shields.io/badge/CRAN-published-brightgreen)](https://CRAN.R-project.org/package=bayesQRsurvey)

Bayesian quantile regression for complex survey data with informative sampling. Supports single-output (MCMC) and multiple-output (EM) estimation with survey weights.

## Installation

```r
# CRAN
install.packages("bayesQRsurvey")

# Development version
devtools::install_github("torodriguezt/bayesQRsurvey")
```

## Usage

### Single-output quantile regression

```r
library(bayesQRsurvey)

fit <- bqr.svy(
  y ~ x1 + x2,
  weights  = data$weight,
  data     = data,
  quantile = 0.5,
  method   = "ald"
)

summary(fit)
plot(fit)
```

### Multiple-output quantile regression

```r
fit_mo <- mo.bqr.svy(
  cbind(y1, y2) ~ x1 + x2,
  weights  = data$weight,
  data     = data,
  quantile = c(0.05, 0.10, 0.25),
  n_dir    = 12
)

summary(fit_mo)
plotQuantileRegion(fit_mo, response = c("y1", "y2"), datafile = data)
```

## References

Nascimento, M. L., & Goncalves, K. C. M. (2024). Bayesian Quantile Regression Models for Complex Survey Data Under Informative Sampling. *Journal of Survey Statistics and Methodology*, 12(4), 1105--1130. [doi:10.1093/jssam/smae015](https://doi.org/10.1093/jssam/smae015)

## Authors

- **Tomas Rodriguez Taborda** -- Universidad Nacional de Colombia (UNAL)
- **Johnatan Cardona Jimenez** -- Universidad Nacional de Colombia (UNAL)
- **Marcus L. Nascimento** -- Getulio Vargas Foundation (FGV EMAp)
- **Kelly Cristina Mota Goncalves** -- Federal University of Rio de Janeiro (UFRJ)

## License

MIT
