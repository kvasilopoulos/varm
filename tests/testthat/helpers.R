dt <- econdata::sw2001[,-1]
obj_dt <- abvar::spec(dt, .endo_lags = 4) %>% varm()
irf_dt <- irf(obj_dt, id = id_chol(), horizon = 24)
