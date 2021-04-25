dt <- econdata::sw2001[,-1]
obj_dt <- varm::spec(dt, .endo_lags = 4) %>% varm()
irf_dt <- irf(obj_dt, id = id_chol(), shock = shock_scale("unit"), horizon = 24)

# autoplot(irf_dt, var_names = c("Inflation", "Unemployment", "Federal Funds Rate"), secondary_bands = sec_bands(shade_color = "red"))

