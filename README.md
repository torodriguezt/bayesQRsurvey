# bayesQRsurvey: Bayesian Weighted Quantile Regression

[![CRAN published](https://img.shields.io/badge/CRAN-published-brightgreen)](https://CRAN.R-project.org/package=bayesQRsurvey)

**bayesQRsurvey** is an R package for **Bayesian weighted quantile regression** for complex survey designs. The package provides both single and multiple-output quantile regression estimation using efficient MCMC and EM algorithms with fast C++ implementations.

---

## ⚙️ Features

- **Survey weights**: Handles complex survey designs with observation weights
- **Multiple algorithms**: MCMC (ALD, Score, Approximate) and EM methods  
- **Fast computation**: C++ implementations using Rcpp, RcppEigen, and RcppArmadillo
- **Comprehensive output**: Detailed summaries with convergence diagnostics
---

## 🚀 Performance Benchmark

The following benchmark compares the R and C++ versions of the main algorithms on a simulated dataset with 100,000 observations and 10 predictors.

| Algorithm            | R Time (seg) | C++ Time (seg) | R RAM  | C++ RAM | Speedup Factor | Memory Saving |
|----------------------|--------------|----------------|--------|---------|----------------|----------------|
| EM_BWQR_AL_MO        | 2.44         | 0.0032         | 2.3 GB | 190 MB  | ×769           | ~12×           |
| MCMC_BWQR_AL         | 12.3         | 0.01           | 2.0 GB | 50 MB   | ×1100          | ~40×           |
| MCMC_BWQR_AP         | 21.78        | 2.81           | ? GB   | ? GB    | ×7.8           | ?              |
| MCMC_BWQR_SL         | 9.7          | 1.4            | ?      | ?       | ×7.1           | ~7.2×          |

> Test environment: R 4.4.2, Windows 11, Intel i5 13600-K, 32 GB RAM

---

## 📦 Installation

### From GitHub (development version)

```r
# Install 'devtools' (or 'remotes') once
install.packages("devtools")

# Install the latest bayesQRsurvey from GitHub
devtools::install_github("torodriguezt/bayesQRsurvey")

# Load the package
library(bayesQRsurvey)
```

---

## 🎯 Main Functions

### Quantile Regression: `bqr.svy()`

Fits Bayesian quantile regression for a single quantile using MCMC methods:

- **ALD (Asymmetric Laplace Distribution)**
- **Score**
- **Approximate**

```r
library(bayesQRsurvey)

# Simulate data
set.seed(123)
n <- 100
x1 <- rnorm(n)
x2 <- runif(n) 
y <- 1 + 2*x1 - 0.5*x2 + rnorm(n)
weights <- runif(n, 0.5, 2)
data <- data.frame(y, x1, x2)

# Fit single quantile regression
fit <- bqr.svy(y ~ x1 + x2, weights = weights, data = data, 
               quantile = 0.5, method = "ald", niter = 5000, burnin = 2500)
summary(fit)
```

### Multiple-output Quantile Regression: `mo.bqr.svy()`

Fits Bayesian quantile regression for multiple quantiles using EM algorithm:

```r
# Fit multiple quantile regression
fit_multi <- mo.bqr.svy(y ~ x1 + x2, weights = weights, data = data, 
                        quantile = c(0.25, 0.5, 0.75), algorithm = 'em')
summary(fit_multi)
```

### Visualization: `plot()`

Create plots showing fitted quantile curves with credible intervals:

```r
# Plot results
plot(fit)
plot(fit_multi)
```

---

## 📚 Methods

The package implements methods based on:

- Marcus L Nascimento, Kelly C M Gonçalves, Bayesian Quantile Regression Models for Complex Survey Data Under Informative Sampling, Journal of Survey Statistics and Methodology, Volume 12, Issue 4, September 2024, Pages 1105–1130, [https://doi.org/10.1093/jssam/smae015](https://academic.oup.com/jssam/article-abstract/12/4/1105/7642687)
- [GitHub Repository: bqr_informative_sampling](https://github.com/marcuslavagnole/bqr_informative_sampling)

---

## 👥 Authors

- **Tomás Rodríguez Taborda**  
  Student, Department of Statistics and Department of Computer and Decision Sciences, Universidad Nacional de Colombia (UNAL)

- **Johnatan Cardona Jiménez**  
  Assistant Professor, Department of Statistics, Universidad Nacional de Colombia (UNAL)
  
- **Marcus L. Nascimento**  
  Postdoctoral Researcher, School of Applied Mathematics, Getulio Vargas Foundation (FGV EMAp)

- **Kelly Cristina Mota Gonçalves**  
  Associate Professor, Department of Statistical Methods, Federal University of Rio de Janeiro (UFRJ)

---

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.
