# Reduced-form OLS, Cholesky IRF (1 std and "1 percent" units), FEVD, forecast, historical decomposition, long-run (Blanchard-Quah)
# varm vs BVAR_ (bvar_.m ols block, iresponse.m, fevd.m, forecasts.m, histdecomp.m, iresponse_longrun.m).  Deterministic comparisons at OLS point estimates.
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm()
set.seed(20261011)
GK <- read_mat5(file.path(EXDIR, "DataGK.mat"))
y <- cbind(logip = GK$logip, logcpi = GK$logcpi, gs1 = GK$gs1, ebp = GK$ebp)       # example_3_irf.m (case 1)
lags <- 12; hor <- 24; K <- 4
cat("data: DataGK.mat [logip logcpi gs1 ebp], T =", nrow(y), " lags =", lags, " hor =", hor, "\n")
# ---- BVAR_ OLS (bvar_.m: varols via rfvar3 with no dummies; Sigma_ols = u'u/(nobs-nk)) ------------
o <- yxb(y, lags); nobs <- nrow(o$Y); nk <- K * lags + 1
Phi <- qr.coef(qr(o$X, tol = 1e-12), o$Y); u <- o$Y - o$X %*% Phi; Sig <- crossprod(u) / (nobs - nk)
fit <- varm(y, p = lags)
rep_("P1 OLS coefficients B equal (tol 1e-6)", mx(fit$B, Phi) < 1e-6, "max|diff|=", signif(mx(fit$B, Phi), 3), " cond(X)=", signif(kappa(o$X, exact = TRUE), 3))
nrm <- function(B) max(abs(crossprod(o$X, o$Y - o$X %*% B)))
rep_("P1b normal-equation residual max|X'(Y-XB)| < 1e-4 (QR/SVD reference gives 1e-7)", nrm(fit$B) < 1e-4, "varm=", signif(nrm(fit$B), 3), " QR=", signif(nrm(Phi), 3), " -> varm OLS solves X'X (cond^2 loss)")
fit0 <- fit; fit$B <- Phi; fit$U <- u; fit$Sigma <- Sig; fit$Sigma_ml <- crossprod(u) / nobs   # engine-level comparisons below use BVAR_'s (accurate) OLS so conditioning does not mask logic
rep_("P2 Sigma (dof-corrected) equal", mx(fit0$Sigma, Sig) < 1e-10, "max|diff|=", signif(mx(fit$Sigma, Sig), 3))
# ---- Cholesky IRF, one std (iresponse.m) ----------------------------------------------------------
ir <- iresponse(Phi[1:(K * lags), ], Sig, hor, diag(K))                      # [var, h, shock]
iv <- irf(fit, id_chol(), horizon = hor - 1, ci = "none")$estimate            # [var, shock, h]
d <- max(sapply(1:hor, function(h) mx(iv[, , h], ir[, h, ])))
rep_("D1 Cholesky IRF (1 std) all vars/shocks/24 horizons", d < 1e-10, "max|diff|=", signif(d, 3))
# ---- "1 percent" IRF (iresponse(...,unit=0): Q*inv(diag(diag(Q)))) vs varm normalise ------------------
irp <- iresponse(Phi[1:(K * lags), ], Sig, hor, diag(K), unit = 0)
dd <- sapply(1:K, function(j) { v <- irf(fit, id_chol(), horizon = hor - 1, ci = "none", normalise = list(shock = j, variable = j, size = 1))$estimate
                                max(sapply(1:hor, function(h) mx(v[, j, h], irp[, h, j]))) })
rep_("D2 unit-impact IRF (BVAR_ unit=0 vs varm normalise)", max(dd) < 1e-10, "max|diff|=", signif(max(dd), 3))
# ---- FEVD (fevd.m: cumulative over hor terms; percent) --------------------------------------------
fe <- fevd_m(hor, Phi[1:(K * lags), ], Sig)
fv <- fevd(fit, id_chol(), horizon = hor - 1, ci = "none")$estimate[, , hor] * 100
rep_("D3 FEVD at horizon 24 (percent; rows var, cols shock)", mx(fe, fv) < 1e-8, "max|diff|=", signif(mx(fe, fv), 3), " row sums:", paste(round(rowSums(fv), 6), collapse = ","))
# ---- forecast (forecasts.m no-shock path; constant only) -------------------------------------------
fh <- 12; init <- y[(nrow(y) - lags + 1):nrow(y), ]
fb <- forecast_noshock(init, Phi, fh, matrix(1, fh, 1)); fvm <- forecast(fit, h = fh, ci = "none")$mean
rep_("F1 point forecast 12 steps", mx(fb, fvm) < 1e-9, "max|diff|=", signif(mx(fb, fvm), 3))
# ---- historical decomposition (histdecomp.m with draw = OLS) -------------------------------------------
hb <- histdecomp(u, Phi, Sig, y, lags); hv <- hd(fit, id_chol())
# BVAR_ innovations whiten with chol(Sigma) = dof-corrected; varm uses the same Sigma for P? compare contributions
dH <- mx(hb$E, hv$contributions)
rep_("H1 HD contributions of each shock (Te x K x K)", dH < 1e-8, "max|diff|=", signif(dH, 3))
rep_("H2 BVAR_ HD adds up to data (its own tol check)", mx(hb$yhat, y[-(1:lags), ]) < 1e-7, "max|diff|=", signif(mx(hb$yhat, y[-(1:lags), ]), 3))
rep_("H3 deterministic+initial-condition part equals varm baseline", mx(hb$B, hv$baseline) < 1e-8, "max|diff|=", signif(mx(hb$B, hv$baseline), 3))
# ---- long run (iresponse_longrun.m) -------------------------------------------------------------------
lr <- iresponse_longrun(Phi[1:(K * lags), ], Sig, hor, lags)
vlr <- irf(fit, id_longrun(), horizon = hor - 1, ci = "none")$estimate
dl <- sapply(1:K, function(j) max(sapply(1:hor, function(h) mx(vlr[, j, h], lr$ir[, h, j]))))
cat("   per-shock max|diff| BVAR_ vs varm long-run IRF:", signif(dl, 3), "\n")
rep_("L1 Blanchard-Quah IRF equal (all 4 shocks)", max(dl) < 1e-8, "max|diff|=", signif(max(dl), 3))

# L2: BVAR_'s `Q(1,1) = Q(1,1)*sign(Q(1,1))` flips only ONE element, not the column.  In ~8% of random variable orderings/datasets
# (35/420 in a scan of Data.mat, DataGK, DataGLP, DataEx) Q(1,1) < 0.  Example: Data.mat [IPI HICP Euribor1Y], p = 1.
Dd <- read_mat5(file.path(EXDIR, "Data.mat")); ye <- cbind(IPI = Dd$IPI, HICP = Dd$HICP, Euribor1Y = Dd$Euribor1Y); pe <- 1; Ke <- 3
oe <- yxb(ye, pe); Pe <- qr.coef(qr(oe$X, tol = 1e-12), oe$Y); ue <- oe$Y - oe$X %*% Pe; Se <- crossprod(ue) / (nrow(oe$Y) - ncol(oe$X))
lre <- iresponse_longrun(Pe[1:(Ke * pe), ], Se, 20, pe)
C1 <- solve(diag(Ke) - t(Pe[1:Ke, ])); Qraw <- solve(C1) %*% t(chol(C1 %*% Se %*% t(C1)))
fe_ <- varm(ye, pe); fe_$B <- Pe; fe_$U <- ue; fe_$Sigma <- Se
vl <- irf(fe_, id_longrun(), horizon = 19, ci = "none")$estimate
cat("   raw Q(1,1) =", signif(Qraw[1, 1], 4), "(negative -> BVAR_ flips only that element)\n")
Qb <- lre$Q
rep_("L2a BVAR_ long-run impact satisfies Q Q' = Sigma when Q(1,1)<0", mx(Qb %*% t(Qb), Se) < 1e-10, "max|QQ'-Sigma|=", signif(mx(Qb %*% t(Qb), Se), 3))
Pv <- vl[, , 1]; rep_("L2b varm long-run impact satisfies P P' = Sigma", mx(Pv %*% t(Pv), Se) < 1e-10, "max|PP'-Sigma|=", signif(mx(Pv %*% t(Pv), Se), 3))
dl2 <- max(sapply(1:20, function(h) mx(vl[, , h], lre$ir[, h, ])))
rep_("L2c long-run IRFs agree up to column sign (Q(1,1)<0 case)", dl2 < 1e-8, "max|diff|=", signif(dl2, 3), "(disagreement is the BVAR_ one-element flip)")
# verify long-run restriction holds in varm for every draw type (C(1)P lower triangular) in the K=4 GK case
Pk <- irf(fit, id_longrun(), horizon = 0, ci = "none")$estimate[, , 1]; Fk <- diag(4) - Reduce(`+`, lapply(1:lags, function(l) t(Phi[(l - 1) * 4 + 1:4, ]))); LRm <- solve(Fk) %*% Pk
rep_("L3 varm: long-run matrix C(1)P is lower triangular", max(abs(LRm[upper.tri(LRm)])) < 1e-9, "max|upper|=", signif(max(abs(LRm[upper.tri(LRm)])), 3))

# ---- H4: historical decomposition with a rotation (histdecomp.m options.Omega) vs varm hd(id_given(P = chol(Sigma) %*% Omega)) ----
Om <- qr.Q(qr(matrix(rnorm(16), 4))); hb2 <- histdecomp(u, Phi, Sig, y, lags, Omega = Om); hv2 <- hd(fit, id_given(t(chol(Sig)) %*% Om))
rep_("H4 HD with an orthogonal rotation of the Cholesky factor", mx(hb2$E, hv2$contributions) < 1e-8, "max|diff|=", signif(mx(hb2$E, hv2$contributions), 3))
# ---- F2: forecast density at a fixed (B, Sigma): covariance of h-step errors = sum_i Theta_i Sigma Theta_i' (forecasts.m uses chol(Sigma) * randn) ----
N <- 100000; Hf <- 6
fd <- varm:::cpp_forecast_draws(array(Phi, c(dim(Phi), N)), array(Sig, c(4, 4, N)), lags, init, matrix(1, Hf, 1), Hf)     # [draw, h, var]
MAr <- varm:::cpp_ma(Phi, 4, lags, Hf - 1)
chk <- sapply(1:Hf, function(h) { Vt <- Reduce(`+`, lapply(1:h, function(i) MAr[, , i] %*% Sig %*% t(MAr[, , i]))); Ve <- cov(fd[, h, ]); max(abs(Ve - Vt) / sqrt(outer(diag(Vt), diag(Vt)))) })
fm <- fb[1:Hf, ]; mz <- max(abs(colMeans(fd[, Hf, ]) - fm[Hf, ]) / sqrt(diag(Sig) * Hf))
rep_("F2a forecast draws: mean = point path (|z|<4.5)", max(abs((apply(fd, c(2, 3), mean) - fm) / sqrt(apply(fd, c(2, 3), var) / N))) < 4.5, "max|z|=", signif(max(abs((apply(fd, c(2, 3), mean) - fm) / sqrt(apply(fd, c(2, 3), var) / N))), 3))
rep_("F2b forecast draws: covariance of h-step errors equals theory (error in correlation units < 1.5%; MC se ~0.5%)", max(chk) < 0.015, "max rel err=", signif(max(chk), 3))
