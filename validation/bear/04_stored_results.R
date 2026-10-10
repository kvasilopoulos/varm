# BEAR-shipped expected results (tests/bear5/results_*.mat, parsed with matread.py) vs varm / vs R replicas of BEAR code.
source("validation/bear/common.R")
set.seed(3)
## ---- A. results_ols.mat: OLS VAR + OLS forecast (data_.xlsx = tutorial data, 1974Q1-2014Q4, 4 lags, const) ----
o <- bear_mat("results_ols.mat", c("data", "Bhat", "sigmahat", "EPS", "forecast_estimates"))
Y <- o$data[, 1:3]; p <- 4; f <- varm(Y, p)
rep_("A1 stored BEAR OLS Bhat vs varm", max(abs(f$B - o$Bhat)) < 1e-10, "maxdiff", format(max(abs(f$B - o$Bhat)), digits = 3))
rep_("A2 stored BEAR OLS sigmahat (T-k) vs varm$Sigma", max(abs(f$Sigma - o$sigmahat)) < 1e-10, "maxdiff", format(max(abs(f$Sigma - o$sigmahat)), digits = 3))
rep_("A3 stored residuals EPS vs varm$U", max(abs(f$U - o$EPS)) < 1e-10)
# BEAR olsforecast.m: point forecast from the sample ending 2013Q4 (Fstartdate=2014q1), 12 periods, with full-sample Bhat
f2 <- f; f2$Y <- Y[1:160, ]
fc <- forecast(f2, h = 12, ci = "none")$mean
fb <- t(sapply(o$forecast_estimates, function(m) m[2, ]))
rep_("A4 OLS point forecast (12 steps) vs BEAR stored", max(abs(t(fc) - fb)) < 1e-9, "maxdiff", format(max(abs(t(fc) - fb)), digits = 3))
# BEAR band: normal, no parameter uncertainty, cumulative sum_{j<=h} Psi_j Sigma Psi_j' (olsforecast.m)
Psi <- bear_irfsim(f$B, diag(3), 3, 1, p, 12)$irf
sf <- Reduce(`+`, lapply(1:12, function(h) Psi[, , h] %*% f$Sigma %*% t(Psi[, , h])), accumulate = TRUE)
lo <- t(sapply(1:12, function(h) fc[h, ] + qnorm(.025) * sqrt(diag(sf[[h]]))))
fl <- t(sapply(o$forecast_estimates, function(m) m[1, ]))
rep_("A5 BEAR forecast band formula (replica) vs stored", max(abs(t(lo) - fl)) < 1e-9)
set.seed(1); fm <- forecast(f2, h = 12, ci = "mc", nboot = 20000, level = .95)
cat("   varm mc 95% lower band vs BEAR normal band (h=1,6,12, var1): varm", round(fm$bands[c(1, 6, 12), 1, 1], 3), " BEAR", round(fl[1, c(1, 6, 12)], 3), "\n")
cat("   (BEAR band = Gaussian, cumulative Psi_j Sigma Psi_j' for j<h: the correct h-step forecast MSE with parameters treated as known)\n")

## ---- B. results_test_data.mat: Minnesota prior 12 (Sigma fixed = diag(sigmahat)), 4 lags, const, ar=.8 l1=.1 l2=.5 l3=1 l4=100 ----
m <- bear_mat("results_test_data.mat", c("data", "Bhat", "arvar", "ar", "beta0", "omega0", "sigma", "betabar", "omegabar", "beta_gibbs", "lambda4"))
Y <- m$data[, 1:3]; mn <- bear_minnesota(Y, 4, TRUE, .8, .1, .5, 1, 100, 12)
cat("   stored lambda4 =", m$lambda4, " ar =", m$ar, "\n")
rep_("B1 BEAR arvar (stored) vs replica of arloop.m", max(abs(mn$arvar - m$arvar)) < 1e-10, "maxdiff", format(max(abs(mn$arvar - m$arvar)), digits = 3))
rep_("B2 BEAR beta0 (stored) vs replica", max(abs(mn$beta0 - m$beta0)) < 1e-12)
rep_("B3 BEAR omega0 (stored) vs replica of mprior.m", max(abs(mn$omega0 - m$omega0)) < 1e-10, "maxdiff", format(max(abs(mn$omega0 - m$omega0)), digits = 3))
rep_("B4 BEAR sigma (stored) vs replica", max(abs(mn$sigma - m$sigma)) < 1e-10)
rep_("B5 BEAR betabar (stored) vs replica of mpost.m", max(abs(mn$betabar - m$betabar)) < 1e-8, "maxdiff", format(max(abs(mn$betabar - m$betabar)), digits = 3))
rep_("B6 BEAR omegabar (stored) vs replica", max(abs(mn$omegabar - m$omegabar)) < 1e-8)
av <- varm:::ar_resid_var(Y, 4)
rep_("B7 varm ar_resid_var vs stored BEAR arvar", max(abs(av / m$arvar - 1)) < 1e-8, "ratio", paste(format(av / m$arvar, digits = 6), collapse = " "), "(df convention)")
cat("   BEAR beta_gibbs mean vs betabar (BEAR self-consistency): max|z|=", round(max(abs(rowMeans(m$beta_gibbs) - m$betabar) / (sqrt(diag(m$omegabar)) / sqrt(ncol(m$beta_gibbs)))), 2), "\n")
## varm has no fixed-Sigma Minnesota prior (prior 11-13): compare nearest equivalents informationally
fo <- varm(Y, 4)
cat("   B8 NOT-TESTED for varm: Minnesota with fixed Sigma (BEAR priors 11/12/13) has no varm counterpart (bvarm is NIW/GLP, bvarm_tri is equation-wise).\n")
cat("   Info: max|betabar_BEAR - OLS| =", format(max(abs(mn$betabar - c(fo$B))), digits = 3), "\n")

## ---- C. CH2019 data: OLS VAR(12) 5 var and GK-style proxy (BEAR irfiv_ols_for_bootstrap_GK.m) vs varm_proxy ----
suppressMessages(library(readxl))
xl <- file.path(BEAR, "BEARX-Toolbox/replications/data/data_CH2019.xlsx")
D <- as.data.frame(suppressMessages(read_excel(xl, "data"))); IVd <- as.data.frame(suppressMessages(read_excel(xl, "IV")))
D <- D[match("1993m1", D$datestr):match("2007m6", D$datestr), ]
ch <- bear_mat("results_test_data_CH2019.mat", c("Bhat", "sigmahat", "Y", "data", "D_estimates"))
Yc <- as.matrix(D[, 2:6]); p <- 12; fc <- varm(Yc, p)
cat("   CH2019 sample rows:", nrow(Yc), " BEAR Y rows:", nrow(ch$Y), "\n")
rep_("C1 CH2019 OLS Bhat (61x5) vs varm", max(abs(fc$B - ch$Bhat)) < 1e-6, "maxdiff", format(max(abs(fc$B - ch$Bhat)), digits = 3))
rep_("C2 CH2019 OLS sigmahat vs varm$Sigma", max(abs(fc$Sigma - ch$sigmahat)) < 1e-8)
z <- IVd$MHF[match(D$datestr, IVd[, 1])]; stopifnot(!anyNA(z))
# BEAR GK: 2SLS impact (u_h on fitted u_1|z), then D1 = b / ||chol(Sigma_IV)^-1 b||, Sigma_IV = EPS'EPS/(T-k)
U <- fc$U; zz <- z[-(1:p)]; k <- 5 * p + 1; Tn <- nrow(U)
fs <- fitted(lm(U[, 1] ~ zz)); b <- c(1, sapply(2:5, function(h) coef(lm(U[, h] ~ fs))[2]))
SigIV <- crossprod(U) / (Tn - k); Cc <- t(chol(SigIV)); qv <- solve(Cc, b); D1 <- b / sqrt(sum(qv^2))
vp <- varm_proxy(fc, m = z, target = 1, size = 1, horizon = 48, nboot = 0)
Th <- bear_irfsim(fc$B, diag(5), 5, 1, p, 49)$irf
ir_b <- t(sapply(1:49, function(h) Th[, , h] %*% (b)))   # IRF with unit impact on var 1 (b[1] = 1)
rep_("C3 proxy IRF, unit impact on target (BEAR GK 2SLS replica vs varm_proxy), 49 horizons", max(abs(vp$estimate - ir_b)) < 1e-4, "maxdiff", format(max(abs(vp$estimate - ir_b)), digits = 3))
cat("   BEAR unit-variance impact D1:", round(D1, 4), "\n   varm B1 (unit-variance):    ", round(drop(vp$B1), 4), "\n   ratio:", round(drop(vp$B1) / D1, 4), "  sqrt((T-k)/T)=", round(sqrt((Tn - k) / Tn), 4), "\n")
rep_("C4 unit-variance B1: ratio varm/BEAR equals sqrt((T-k)/T) exactly (convention: BEAR Sigma/(T-k), varm Sigma/T)", max(abs(drop(vp$B1) / D1 - sqrt((Tn - k) / Tn))) < 1e-4)
cat("   BEAR-stored CH2019 posterior-median D[,1] (Bayesian proxy, IRFt=5):", round(ch$D_estimates[1:5], 3), "\n   frequentist proxy B1 (varm) sign/scale for reference:", round(drop(vp$B1), 3), "\n")
cat("   C5 NOT-TESTED: BEAR IRFt=5 is the Caldara-Herbst (2019) Metropolis-within-Gibbs Bayesian proxy-SVAR; varm documents that it is NOT implemented (varm_proxy_bayes/arw/gkr are different posteriors).\n")
