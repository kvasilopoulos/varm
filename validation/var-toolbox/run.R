# Compare varm with R ports of Cesa-Bianchi's VAR-Toolbox. Run from the varm root:
#   Rscript validation/var-toolbox/run.R            (about 2-4 minutes)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("validation/var-toolbox/mlref.R"); source("validation/var-toolbox/load.R")
RES <- data.frame(test = character(), stat = character(), value = numeric(), tol = numeric(), status = character())
rec <- function(test, stat, value, tol, note = "") {
  st <- if (is.na(value)) "NA" else if (value <= tol) "PASS" else "FAIL"
  RES[nrow(RES) + 1, ] <<- list(test, stat, value, tol, st)
  cat(sprintf("%-34s %-30s %12.4g (tol %.0e) %s %s\n", test, stat, value, tol, st, note))
}
mx <- function(a, b) max(abs(a - b))
# varm structural array [resp, shock, h+1] -> MATLAB layout [h, resp, shock]
tml <- function(a) aperm(unname(a), c(3, 1, 2))
bq <- rd("BQ1989"); uh <- rd("Uhlig2005"); gk <- rd("GK2015"); gk[gk == 123456789] <- NA
gkY <- gk[, c("gs1", "logcpi", "logip", "ebp")]; gkZ <- gk[, "ff4_tc"]

## ---- A. reduced-form estimation ------------------------------------------------
cat("\n== A. estimation ==\n")
for (nm in c("BQ(p=8)", "Uhlig(p=12)", "GK(p=12)")) {
  D <- switch(nm, "BQ(p=8)" = bq, "Uhlig(p=12)" = uh, "GK(p=12)" = gkY); p <- if (nm == "BQ(p=8)") 8 else 12
  f <- varm(D, p); v <- ml_var(D, p, 1); K <- ncol(D)
  Bm <- v$Ft[c((1 + 1):(K * p + 1), 1), ]            # MATLAB [const, lags] -> varm [lags, const]
  qrB <- qr.solve(v$X, v$Y)[c(2:(K * p + 1), 1), ]   # QR solution = numerically safest reference
  rec(nm, sprintf("coef vs ML normal-eq (cond=1e%.0f)", log10(kappa(crossprod(v$X), exact = TRUE))), mx(unname(f$B), Bm), 1e-5)
  rec(nm, "coef vs QR solve", mx(unname(f$B), qrB), if (nm == "BQ(p=8)") 1e-8 else 1e-5)
  rec(nm, "ML normal-eq vs QR (ML's own error)", mx(Bm, qrB), 1e-5)
  rec(nm, "Sigma (dof-corrected)", mx(unname(f$Sigma), v$sigma), 1e-10)
  rec(nm, "residuals", mx(unname(f$U), v$resid), 1e-5)
  rec(nm, "max eigenvalue", abs(max(stability(f)) - v$maxEig), 1e-6)
}
fq <- varm(bq, 4, trend = TRUE); vq <- ml_var(bq, 4, 2); K <- 2
rec("BQ const+trend p=4", "lag coef", mx(unname(fq$B[1:8, ]), vq$Ft[3:10, ]), 1e-8)
rec("BQ const+trend p=4", "trend slope", mx(unname(fq$B["trend", ]), vq$Ft[2, ]), 1e-8)
rec("BQ const+trend p=4", "intercept (varm trend=1..T vs ML 1..Te)", mx(unname(fq$B["const", ]), vq$Ft[1, ]), 1e-8,
    "<- expected shift = p*slope")
rec("BQ const+trend p=4", "c_varm - c_ML + p*slope (should be 0)", mx(unname(fq$B["const", ] - vq$Ft[1, ]), -4 * vq$Ft[2, ]), 1e-8)
rec("BQ const+trend p=4", "residuals", mx(unname(fq$U), vq$resid), 1e-8)

## ---- B. Cholesky IRF -------------------------------------------------------------
cat("\n== B. Cholesky IRF ==\n")
H <- 48
f <- varm(gkY, 12); v <- ml_var(gkY, 12); B <- ml_B(v, "short")
ir <- irf(f, id_chol(), horizon = H, ci = "none")$estimate
rec("GK chol IRF 1sd", "max|diff| vs ml (wold)", mx(tml(ir), ml_ir(v, B, H + 1)), 1e-6)
rec("GK chol IRF 1sd", "max|diff| vs ml (companion)", mx(tml(ir), ml_ir(v, B, H + 1, recurs = "comp")), 1e-6)
ir1 <- irf(f, id_chol(), horizon = H, ci = "none", normalise = list(shock = "gs1", variable = "gs1", size = 1))$estimate
ml1 <- ml_ir(v, B, H + 1, impact = 1)
rec("GK chol IRF unit (impact=1)", "shock gs1 (normalise)", mx(tml(ir1)[, , 1], ml1[, , 1]), 1e-6)
# ordering: MATLAB users reorder columns of data; varm has order=
ord <- c(3, 1, 4, 2); fo <- varm(gkY, 12); iro <- irf(fo, id_chol(order = ord), horizon = 8, ci = "none")$estimate
vo <- ml_var(gkY[, ord], 12); mo <- ml_ir(vo, ml_B(vo, "short"), 9)
rec("GK chol, order=c(3,1,4,2)", "vs ML on permuted data", mx(tml(iro)[, ord, ord], mo), 1e-8)
f2 <- varm(uh, 12); v2 <- ml_var(uh, 12)
rec("Uhlig chol IRF (6 var,p=12)", "max|diff|", mx(tml(irf(f2, horizon = 60, ci = "none")$estimate), ml_ir(v2, ml_B(v2, "short"), 61)), 1e-5)

## ---- C. long-run (Blanchard-Quah) --------------------------------------------------
cat("\n== C. long-run ==\n")
fb <- varm(bq, 8); vb <- ml_var(bq, 8); Bl <- ml_B(vb, "long")
irl <- irf(fb, id_longrun(), horizon = 39, ci = "none")$estimate
rec("BQ long-run IRF", "max|diff| vs ml", mx(tml(irl), ml_ir(vb, Bl, 40)), 1e-8)
rec("BQ long-run IRF", "impact matrix", mx(unname(irl[, , 1]), Bl), 1e-8)
A1 <- diag(2) - Reduce(`+`, lapply(1:8, function(l) t(fb$B[(l - 1) * 2 + 1:2, ])))
LRm <- solve(A1) %*% irl[, , 1]
rec("BQ long-run", "LR matrix upper elem (should be 0)", abs(LRm[1, 2]), 1e-10)
rec("BQ long-run", "B B' - Sigma", mx(irl[, , 1] %*% t(irl[, , 1]), fb$Sigma), 1e-10)

## ---- D. FEVD -----------------------------------------------------------------------
cat("\n== D. FEVD ==\n")
fv <- fevd(f, id_chol(), horizon = H, ci = "none")$estimate         # [resp, shock, h]
vd <- ml_vd(v, B, H + 1)                                            # [h, shock, resp] percent
rec("GK chol FEVD", "max|diff| (percentage points)", mx(100 * aperm(unname(fv), c(3, 2, 1)), vd), 1e-3)
fvl <- fevd(fb, id_longrun(), horizon = 39, ci = "none")$estimate
rec("BQ long-run FEVD", "max|diff| (percentage points)", mx(100 * aperm(unname(fvl), c(3, 2, 1)), ml_vd(vb, Bl, 40)), 1e-6)
rec("BQ long-run FEVD", "sum over shocks == 1", mx(apply(fvl, c(1, 3), sum), 1), 1e-10)

## ---- E. historical decomposition ----------------------------------------------------
cat("\n== E. HD ==\n")
hdv <- hd(f, id_chol()); hm <- ml_hd(v, B); p <- 12
sh <- hm$shock[-(1:p), , ]                                          # [t, shock, var]
rec("GK chol HD shocks", "contributions max|diff|", mx(aperm(hdv$contributions, c(1, 3, 2)), sh), 1e-3)
rec("GK chol HD", "baseline vs ML init+const", mx(hdv$baseline, hm$init + hm$const), 1e-3)
hl <- hd(fb, id_longrun()); vm <- ml_hd(vb, Bl)
rec("BQ long-run HD shocks", "contributions max|diff|", mx(aperm(hl$contributions, c(1, 3, 2)), vm$shock[-(1:8), , ]), 1e-6)

## ---- F. proxy / IV (Gertler-Karadi data) ----------------------------------------------
cat("\n== F. IV (GK) ==\n")
Bm <- ml_iv(v, gkZ)                                                  # unit-sd impact, gs1 first
Zf <- gkZ; Zf[is.na(Zf)] <- 0                                        # varm needs zero-fill
vp <- varm_proxy(f, Zf, target = "gs1", size = 1, horizon = H, nboot = 0)
mlIR <- ml_ir(v, cbind(Bm, matrix(0, 4, 3)), H + 1)[, , 1]           # [h, var]
irN <- t(sapply(1:(H + 1), function(i) mlIR[i, ] / Bm[1]))           # normalise gs1 impact = 1
rec("GK IV (zero-fill) vs ML NaN-drop", "normalised IRF max|diff|", mx(unname(vp$estimate), irN), 1e-8,
    "(ML drops NaN rows, centres; varm uses uncentred sum over full sample)")
rec("GK IV", "gamma-direction cosine diff", 1 - abs(sum(vp$gamma * Bm) / sqrt(sum(vp$gamma^2) * sum(Bm^2))), 1e-3)
rec("GK IV", "B1 (unit-var) vs ML Biv rel.diff", mx(vp$B1[, 1], Bm) / max(abs(Bm)), 5e-2,
    sprintf("varm=%s ML=%s", paste(round(vp$B1[, 1], 4), collapse = ","), paste(round(Bm, 4), collapse = ",")))
# demean the instrument on its observed sample (removes the centring difference)
Zc <- gkZ; Zc[!is.na(Zc)] <- Zc[!is.na(Zc)] - mean(Zc, na.rm = TRUE); Zc[is.na(Zc)] <- 0
vpc <- varm_proxy(f, Zc, target = "gs1", size = 1, horizon = H, nboot = 0)
rec("GK IV, instrument demeaned", "normalised IRF max|diff|", mx(unname(vpc$estimate), irN), 1e-6)
idp <- irf(f, id_proxy(Zf, f, target = "gs1"), horizon = H, ci = "none")$estimate[, 1, ]
rec("GK id_proxy vs varm_proxy", "max|diff|", mx(t(unname(idp)), unname(vp$estimate)), 1e-10)
# FEVD share of identified shock: ML uses total variance denominators
vdm <- ml_vd(v, cbind(Bm, matrix(0, 4, 3)), H + 1)[, 1, ]
fvp <- fevd(f, id_proxy(Zf, f, target = "gs1", size = Bm[1]), horizon = H, ci = "none")$estimate[, 1, ]
rec("GK IV FEVD (h=48, gs1 from MP shock)", "ML % vs varm %", abs(vdm[H + 1, 1] - 100 * fvp[1, H + 1]), 1,
    sprintf("ML=%.2f%% varm=%.2f%%", vdm[H + 1, 1], 100 * fvp[1, H + 1]))
# consistent FEVD: varm proxy IRF (size=ML unit-var impact) over total variance from the Cholesky-equivalent MSE
irP <- irf(f, id_proxy(Zf, f, target = "gs1", size = Bm[1]), horizon = H, ci = "none")$estimate[, 1, ]
irC <- irf(f, id_chol(), horizon = H, ci = "none")$estimate
tot <- apply(irC^2, c(1, 3), sum); tot <- t(apply(tot, 1, cumsum))
shr <- 100 * t(apply(irP^2, 1, cumsum)) / tot
rec("GK IV FEVD recomputed w/ total variance", "varm-IRF based vs ML (pp, h=48, 4 vars)", mx(shr[, H + 1], vdm[H + 1, ]), 6,
    sprintf("diff only from B1 scale/centring; ML=%s", paste(round(vdm[H + 1, ], 2), collapse = ",")))
rel <- max(abs(vp$estimate - irN)) / max(abs(irN)); cat(sprintf("  normalised-IRF diff as share of max|IRF|: %.4f\n", rel))
# sign of impact when first-stage slope negative
vn <- varm_proxy(f, -Zf, target = "gs1", size = 1, horizon = 4, nboot = 0)
rec("GK IV, flipped instrument", "normalised IRF unchanged", mx(unname(vn$estimate), unname(vp$estimate[1:5, ])), 1e-10)

## ---- G. bootstrap bands -----------------------------------------------------------------
cat("\n== G. bootstrap (BQ, long-run, p=8, 3000 draws each) ==\n")
set.seed(11); nb <- 3000
mlb <- ml_boot(vb, "long", 40, nb, "bs")
dv <- draws(fb, "bootstrap", nb); sv <- varm:::structural_irf(dv$B, dv$Sigma, 2, 8, 39L, id_longrun(), c("y", "u"))$irf
vv <- aperm(sv, c(3, 1, 2, 4))                                       # [h, resp, shock, draw]
q95 <- function(a) apply(a, 1:3, ml_prctile, c(.025, .975))
qm <- q95(mlb); qv <- q95(vv)
w <- function(q) q[2, , , ] - q[1, , , ]
wm <- w(qm); wv <- w(qv)
rec("BQ long-run bootstrap 95%", "mean band-width ratio varm/ML", abs(mean(wv / wm, na.rm = TRUE) - 1), 0.02,
    sprintf("ratio=%.3f", mean(wv / wm, na.rm = TRUE)))
rec("BQ long-run bootstrap", "draw-mean IRF max|diff| (in sd units)",
    max(abs(apply(vv, 1:3, mean) - apply(mlb, 1:3, mean)) / apply(mlb, 1:3, sd)), 0.2)
pm <- which(abs(apply(vv, 1:3, mean) - apply(mlb, 1:3, mean)) / apply(mlb, 1:3, sd) == max(abs(apply(vv, 1:3, mean) - apply(mlb, 1:3, mean)) / apply(mlb, 1:3, sd)), arr.ind = TRUE)
cat("  worst mean-diff cell [h,resp,shock]:", pm[1, ], "\n")
mlb2 <- ml_boot(vb, "long", 40, nb, "bs", scale = TRUE)      # same algorithm plus varm's sqrt(Te/(Te-m)) rescale
rec("BQ long-run bootstrap + ML port with varm's dof rescale", "band-width ratio varm/ML", abs(mean(wv / w(q95(mlb2)), na.rm = TRUE) - 1), 0.02,
    sprintf("ratio=%.3f", mean(wv / w(q95(mlb2)), na.rm = TRUE)))
rec("BQ long-run bootstrap + dof rescale", "draw-mean IRF max|diff| (sd units)",
    max(abs(apply(vv, 1:3, mean) - apply(mlb2, 1:3, mean)) / apply(mlb2, 1:3, sd)), 0.2)
cat("  dof-rescale factor sqrt(Te/(Te-m)) used by varm bootstrap:", sqrt(fb$nobs / (fb$nobs - fb$npar)), "\n")
# faithful MATLAB behaviour inside varm? simulate varm with unscaled residuals is not possible w/o editing R/, report ratio only
set.seed(12); mlw <- ml_boot(ml_var(gkY, 12), "short", 13, 600, "wild")
dvw <- draws(f, "wild", 600); svw <- aperm(varm:::structural_irf(dvw$B, dvw$Sigma, 4, 12, 12L, id_chol(), colnames(gkY))$irf, c(3, 1, 2, 4))
rec("GK chol wild bootstrap", "band-width ratio varm/ML", abs(mean(w(q95(svw)) / w(q95(mlw)), na.rm = TRUE) - 1), 0.05,
    sprintf("ratio=%.3f", mean(w(q95(svw)) / w(q95(mlw)), na.rm = TRUE)))
# quantile type
t5 <- apply(vv, 1:3, quantile, c(.025, .975), type = 5); t7 <- apply(vv, 1:3, quantile, c(.025, .975), type = 7)
rec("quantile rule (N=3000)", "type 5 (MATLAB) vs type 7 (varm) band end, sd units", max(abs(t5 - t7) / rep(apply(vv, 1:3, sd), each = 2)), 0.05, "convention")
# is explosive-draw screening implemented?
dx <- draws(f, "bootstrap", 1000); ex <- mean(apply(stability(dx), 1, max) >= .9999)
cat(sprintf("  explosive/unit-root bootstrap draws kept by varm (GK p=12): %.1f%% (MATLAB would reject them)\n", 100 * ex))
dxb <- draws(fb, "bootstrap", 1000); cat(sprintf("  same for BQ p=8: %.1f%%\n", 100 * mean(apply(stability(dxb), 1, max) >= .9999)))

## ---- H. sign restrictions (Uhlig 2005 replication design) ---------------------------------
cat("\n== H. sign restrictions (Uhlig data, p=12, 6 vars, sr_hor=6, point estimate) ==\n")
SIGN <- matrix(0, 6, 6); SIGN[2, 1] <- -1; SIGN[3, 1] <- -1; SIGN[5, 1] <- -1; SIGN[6, 1] <- 1
set.seed(21); nacc <- 3000; Bs <- list(); tried <- 0
while (length(Bs) < nacc) { o <- ml_sign_once(v2, SIGN, 6, 500); tried <- tried + o$tried; if (!is.null(o$B)) Bs[[length(Bs) + 1]] <- o$B }
IRm <- sapply(Bs, function(B) ml_ir(v2, B, 61)[, , 1], simplify = "array")      # [h, var, draw]
rownames(SIGN) <- colnames(uh); colnames(SIGN) <- paste0("s", 1:6)
set.seed(22); idS <- id_sign(SIGN, horizon = 6, nrot = nacc, maxtry = 200000L)
ri <- irf(f2, idS, horizon = 60, ci = "none", level = 0.68)
rec("Uhlig sign: acceptance", "ML accept/rot vs varm accept", abs(nacc / tried - ri$accept_rate), 1,
    sprintf("ML=%.4f varm=%.4f (ML scans all columns, varm tests column 1 only)", nacc / tried, ri$accept_rate))
IRv <- ri$draws[, 1, , ]                                              # [var, h, draw]
mm <- apply(IRm, 1:2, median); mv <- t(apply(IRv, 1:2, median))
sdIR <- apply(IRm, 1:2, sd)
rec("Uhlig sign: median IRF", "max|diff| in rotation-sd units", max(abs(mm - mv) / sdIR), 0.15)
b <- function(a, pr) apply(a, 1:2, ml_prctile, pr)
rec("Uhlig sign: 16/84% bands", "max|diff| in rotation-sd units",
    max(abs(b(IRm, .16) - t(apply(IRv, 1:2, ml_prctile, .16))) / sdIR, abs(b(IRm, .84) - t(apply(IRv, 1:2, ml_prctile, .84))) / sdIR), 0.2)
mn <- function(a) apply(a, 1:2, min); mxv <- function(a) apply(a, 1:2, max)
rec("Uhlig sign: 2.5/97.5% of rotations", "max|diff| in rotation-sd units",
    max(abs(b(IRm, .025) - t(apply(IRv, 1:2, ml_prctile, .025))), abs(b(IRm, .975) - t(apply(IRv, 1:2, ml_prctile, .975)))) / mean(sdIR), 0.4)
viol <- sum(sapply(1:200, function(i) any(ri$draws[c(2, 3, 5), 1, 1:6, i] > 0) || any(ri$draws[6, 1, 1:6, i] < 0)))
rec("Uhlig sign: restrictions hold in draws", "violations (of 200)", viol, 0)

## ---- I. posterior draw dispersion (VARdrawpost vs draws(mc)) ----------------------------------
cat("\n== I. parameter-uncertainty draws ==\n")
set.seed(31); nd <- 4000; Te <- fb$nobs; Si <- solve(fb$Sigma)
wr <- function(S, df) { L <- t(chol(S)); A <- matrix(0, 2, 2); for (i in 1:2) { A[i, i] <- sqrt(rchisq(1, df - i + 1)); if (i > 1) for (j in 1:(i - 1)) A[i, j] <- rnorm(1) }; LA <- L %*% A; LA %*% t(LA) }
SigM <- sapply(1:nd, function(i) { s <- solve(wr(Si / Te, Te)); s[1, 1] })     # MATLAB VARdrawpost sigma_draw[1,1]
dm <- draws(fb, "mc", nd)
rec("BQ draws: Sigma[1,1] mean", "MATLAB-style vs varm mc rel.diff", abs(mean(SigM) - mean(dm$Sigma[1, 1, ])) / fb$Sigma[1, 1], 0.03)
rec("BQ draws: Sigma[1,1] sd", "ratio varm/ML", abs(sd(dm$Sigma[1, 1, ]) / sd(SigM) - 1), 0.08,
    sprintf("ratio=%.3f (ML df=T, varm df=T-m; theory %.3f)", sd(dm$Sigma[1, 1, ]) / sd(SigM), sqrt(Te / (Te - fb$npar))))

## ---- J. local projections -----------------------------------------------------------------------
cat("\n== J. LP ==\n")
Pm <- rd("Primer"); Pd <- cbind(dlngdp = 100 * c(NA, diff(log(Pm[, "gdp"]))), i1yr = c(NA, diff(Pm[, "i1yr"]))); Pd <- Pd[-1, ]
Hh <- 12; nl <- 4
ml <- ml_lp(Pd[, 1], Pd[, 2], Pd, nl, Hh + 1)                          # MATLAB LPmodel, impact=1, h=1..H+1 -> varm h=0..H
vl <- lp(Pd, shock = "i1yr", horizon = Hh, lags = nl, nw_lag = NULL)
rec("LP-OLS (Primer)", "beta max|diff|", mx(vl$estimate[, "dlngdp"], ml[, 1]), 1e-8)
rec("LP-OLS (Primer)", "NW se: varm default (L=h+1) vs ML (L=h)", mx(vl$se[, "dlngdp"], ml[, 2]), 1e-8, "(bandwidth convention)")
nwl <- sapply(0:Hh, function(h) lp(Pd, shock = "i1yr", horizon = Hh, lags = nl, nw_lag = h)$se[h + 1, "dlngdp"])
rec("LP-OLS (Primer)", "NW se with nw_lag=h", mx(nwl, ml[, 2]), 1e-8)
# LP-IV on GK data (rows with observed ff4)
ok <- which(!is.na(gkZ)); G2 <- cbind(gkY[ok, ], z = gkZ[ok]); G2d <- G2[, c("gs1", "logip", "ebp")]
mli <- ml_lp(G2d[, "logip"], G2d[, "gs1"], G2d, 6, 13, iv = G2[, "z"])
nwi <- t(sapply(0:12, function(h) { r <- lp(G2d, shock = "gs1", horizon = 12, lags = 6, iv = G2[, "z"], nw_lag = h); c(r$estimate[h + 1, "logip"], r$se[h + 1, "logip"]) }))
rec("LP-IV (GK, FF4 -> gs1)", "beta max|diff|", mx(nwi[, 1], mli[, 1]), 1e-8)
rec("LP-IV (GK, FF4 -> gs1)", "se max|diff| (nw_lag=h)", mx(nwi[, 2], mli[, 2]), 1e-8)

write.csv(RES, "validation/var-toolbox/results.csv", row.names = FALSE)
cat("\nSummary:", table(RES$status), "\n"); print(table(RES$status))
