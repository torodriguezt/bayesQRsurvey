devtools::load_all()

plotQuantileRegion(fit_mo3, response = c("wgt", "hgt"),
                   datafile = Anthro, xValue = c(1, 20, 20^2, 1))

print(fit_mo3)
summary(fit_mo3)