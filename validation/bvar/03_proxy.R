# Proxy-SVAR (external instrument): varm vs BVAR_ iresponse_proxy.m (Mertens-Ravn/Gertler-Karadi 2SLS), Gertler-Karadi data of example_3_irf.m (case 6)
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm(); suppressMessages(library(readxl))
set.seed(20261011)
GK <- read_mat5(file.path(EXDIR, "DataGK.mat"))
y <- cbind(gs1 = GK$gs1, logip = GK$logip, logcpi = GK$logcpi, ebp = GK$ebp); lags <- 12; hor <- 24; K <- 4
fd <- as.data.frame(read_excel(file.path(EXDIR, "factor_data.xlsx"), sheet = "factor_data")); prox <- fd$ff4_tc   # numi(:,4)
Tm <- length(prox); T <- nrow(y)
cat("data: DataGK [gs1 logip logcpi ebp], T =", T, "; proxy ff4_tc (factor_data.xlsx col 4), T_m =", Tm, " aligned to the LAST", Tm, "residuals (T_m_end = 0)\n")
o <- yxb(y, lags); nk <- K * lags + 1; Phi <- qr.coef(qr(o$X, tol = 1e-12), o$Y); u <- o$Y - o$X %*% Phi; Sig <- crossprod(u) / (nrow(o$Y) - nk)
fit <- varm(y, lags); fit$B <- Phi; fit$U <- u; fit$Sigma <- Sig       # accurate OLS (varm's normal-equation OLS differs at 6e-8, see 02_dynamics.R)
bm <- iresponse_proxy(u, prox, Phi, Sig, lags, hor)                     # BVAR_ replica
cat("   BVAR_ b1 (unit-variance impact vector):", signif(bm$b1, 5), "\n")
# ---- X1: the first-stage / 2SLS direction: varm gamma-direction equals BVAR_ b21/b11 only if the proxy is centred over its sample ----------
mfull <- c(rep(0, T - Tm), prox)                                         # zero-filled, as varm requires
vp <- varm_proxy(fit, mfull, target = 1, size = bm$b1[1], horizon = hor - 1, nboot = 0)
d_raw <- max(abs(vp$estimate - bm$irs)); rel_raw <- max(abs(vp$estimate[1, ] / vp$estimate[1, 1] - bm$b1 / bm$b1[1]))
rep_("X1 varm_proxy (zero-filled, uncentred) vs BVAR_ (sample-centred) IRF", d_raw < 1e-6, "max|diff|=", signif(d_raw, 3), "  impact-ratio diff=", signif(rel_raw, 3), " (centring of residuals/proxy over the proxy sample)")
mc <- c(rep(0, T - Tm), prox - mean(prox))                                # centre the proxy; residual centring needs a zero-mean u in the sample
# BVAR_ regresses with an intercept in the proxy sample => uses Cov(u - mean_s(u), m - mean(m)); with m centred and zero outside the sample, sum(u*m) = sum((u - ubar_s) * (m - mbar)) exactly.
vc <- varm_proxy(fit, mc, target = 1, size = bm$b1[1], horizon = hor - 1, nboot = 0)
d_c <- max(abs(vc$estimate - bm$irs))
rep_("X2 same, with the proxy centred on its sample (exact match expected)", d_c < 1e-8, "max|diff|=", signif(d_c, 3))
# ---- X3: unit-variance scale b11 (BVAR_ uses sqrt(Sigma11 - b12b12')): varm proxy_B1 formula ---------------------------------------------
ix <- (T - lags - Tm + 1):(T - lags); Sum_c <- crossprod(u[ix, ], prox - mean(prox)) / Tm; Smm_c <- sum((prox - mean(prox))^2) / Tm
B1v <- varm:::proxy_B1(Sig, Sum_c, matrix(Smm_c), 1L)
rep_("X3 varm:::proxy_B1 (Mertens-Ravn unit-variance impact) vs BVAR_ b1", mx(drop(B1v), bm$b1) < 1e-8, "max|diff|=", signif(mx(drop(B1v), bm$b1), 3))
# ---- X4: the published varm_proxy B1 uses ML Sigma and uncentred moments: size of the convention gap ---------------------------------
cat("   varm_proxy$B1 (ML Sigma, uncentred, zero-filled):", signif(drop(vp$B1), 5), "\n")
cat("   relative gap in b11:", signif(abs(vp$B1[1] / bm$b1[1] - 1), 3), " (dof n/(n-k) in Sigma and the centring convention)\n")
# ---- X5: id_proxy(): same direction/normalisation ----------------------------------------------------------------------------------------------
ip <- id_proxy(mfull, fit, target = 1, size = bm$b1[1]); ri <- irf(fit, ip, horizon = hor - 1, ci = "none")$estimate[, 1, ]
rep_("X5 irf(id_proxy) == varm_proxy (r = 1)", mx(t(ri), vp$estimate) < 1e-10, "max|diff|=", signif(mx(t(ri), vp$estimate), 3))
# ---- X6: first-stage statistics on a sample where the proxy covers every residual (no zero-fill) -----------------------------------------------
ys <- y[(T - Tm - lags + 1):T, ]; os <- yxb(ys, lags); Ps <- qr.coef(qr(os$X, tol = 1e-12), os$Y); us <- os$Y - os$X %*% Ps; Ss <- crossprod(us) / (nrow(os$Y) - nk)
fs <- varm(ys, lags); fs$B <- Ps; fs$U <- us; fs$Sigma <- Ss
ms <- c(rep(0, lags), prox)
Z <- cbind(1, prox); Pb <- solve(crossprod(Z), crossprod(Z, us)); Rm <- us - Z %*% Pb; Rc <- us[, 1] - mean(us[, 1])
SST <- sum(Rc^2); SSE <- sum(Rm[, 1]^2); Fm <- ((SST - SSE) / 1) / (SSE / (Tm - 2)); R2 <- 1 - SSE / SST
S0 <- crossprod(Z * Rm[, 1]) / Tm; Av <- solve(crossprod(Z) / Tm) %*% S0 %*% solve(crossprod(Z) / Tm); Wm <- Tm * Pb[2, 1]^2 / Av[2, 2]
fv <- proxy_fstat(fs, ms, 1)
rep_("X6a first-stage F (BVAR_ F_m vs varm F)", abs(fv[["F"]] - Fm) / Fm < 1e-8, "BVAR_=", signif(Fm, 6), " varm=", signif(fv[["F"]], 6))
cat("   robust F: BVAR_ (White HC0 sandwich) =", signif(Wm, 6), " varm F_robust (sandwich HC1) =", signif(fv[["F_robust"]], 6), " ratio=", signif(fv[["F_robust"]] / Wm, 6), " (HC1/HC0 = (T-2)/T =", signif((Tm - 2) / Tm, 6), "inverse)\n")
rep_("X6b robust F equal up to HC1 vs HC0 small-sample factor", abs(fv[["F_robust"]] * ((Tm - 2) / Tm)^-1 / Wm - 1) < 1e-6 || abs(fv[["F_robust"]] / Wm - Tm / (Tm - 2)) < 1e-6 || abs(fv[["F_robust"]] / Wm - (Tm - 2) / Tm) < 1e-6, "")
# ---- X7: same on the subsample (proxy covers all residuals): IRFs exactly equal even without centring ------------------------------------------
bs <- iresponse_proxy(us, prox, Ps, Ss, lags, hor); vs <- varm_proxy(fs, ms, target = 1, size = bs$b1[1], horizon = hor - 1, nboot = 0)
rep_("X7 sub-sample with full proxy coverage: IRF direction equal (VAR residual means are 0 so centring is immaterial for u; proxy still uncentred)", max(abs(vs$estimate - bs$irs)) < 1e-6, "max|diff|=", signif(max(abs(vs$estimate - bs$irs)), 3), "  impact-ratio diff =", signif(max(abs(vs$estimate[1, ] / vs$estimate[1, 1] - bs$b1 / bs$b1[1])), 3))
