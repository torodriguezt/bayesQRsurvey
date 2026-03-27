library(bayesQRsurvey)
library(ggplot2)
data("Anthro")
str(Anthro)
Anthro = Anthro[-which(Anthro$wgt>30),]
Anthro$age2 = Anthro$age**2
Anthro$sex2 = ifelse(Anthro$sex==1, 1, 0)
Anthro$SEX = as.factor(ifelse(Anthro$sex==1, 1, 0))
Anthro = na.omit(Anthro)

set.seed(50) 

fit_ald <- bqr.svy(wgt ~ age + age2 + sex, weights = Anthro$dweight, data = Anthro, 
                   niter = 60000, thin = 1, burnin = 30000)

summary(fit_ald, digits = 3)

fit_score <- bqr.svy(wgt ~ age + age2 + sex2, weights = Anthro$dweight, 
                     data = Anthro, method = "score", niter = 300000, thin = 2, burnin = 200000)


summary(fit_score, digits = 3)


Anthro$AGE =  Anthro$age - mean(Anthro$age)
Anthro$AGE2= Anthro$AGE**2
set.seed(50) #Set the seed to get the same results (optional)
fit_score <- bqr.svy(wgt ~ AGE + AGE2 + sex2, weights = Anthro$dweights, data = Anthro,
                     method = "score", niter = 300000, thin = 2, burnin = 200000)

summary(fit_score, digits = 3)

# =====================================================
#  Multiple-Output Bayesian Quantile Regression (mo.bqr.svy)
#  Bivariate response: wgt and hgt
# =====================================================

# Directions: 2D unit vectors (d=2, K=4 directions)
U <- matrix(c(1, 0,
              0, 1,
              1/sqrt(2), 1/sqrt(2),
              1/sqrt(2), -1/sqrt(2)),
            nrow = 2, ncol = 4)

# Orthogonal complements: each is a 2x1 matrix orthogonal to the corresponding column of U
gamma_U <- list(
  matrix(c(0, 1), ncol = 1),
  matrix(c(1, 0), ncol = 1),
  matrix(c(1/sqrt(2), -1/sqrt(2)), ncol = 1),
  matrix(c(1/sqrt(2),  1/sqrt(2)), ncol = 1)
)

set.seed(50)
fit_mo <- mo.bqr.svy(
  cbind(wgt, hgt) ~ age + age2 + sex2,
  weights = dweight,
  data = Anthro,
  quantile = 0.5,
  U = U,
  gamma_U = gamma_U
)

fit_mo
summary(fit_mo)

# Using n_dir instead of specifying U and gamma_U manually
set.seed(50)
fit_mo2 <- mo.bqr.svy(
  cbind(wgt, hgt) ~ age + age2 + sex2,
  weights = dweight,
  data = Anthro,
  quantile = 0.5,
  n_dir = 4
)

fit_mo2
summary(fit_mo2)



444444444