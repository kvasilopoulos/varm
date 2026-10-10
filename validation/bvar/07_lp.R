# Local projections: varm lp() vs BVAR_ directmethods.m / ols_reg.m (robust_se_ = 1, Hamilton Newey-West with L = lags + h + 1)
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm()
set.seed(20261011)
GK <- read_mat5(file.path(EXDIR, "DataGK.mat")); y <- cbind(logip = GK$logip, logcpi = GK$logcpi, gs1 = GK$gs1, ebp = GK$ebp); K <- 4; lags <- 4; H <- 12; T <- nrow(y)
# ols_reg.m with robust_se_ = 1 (Hamilton 10.5.20): returns beta (K x ny) and se
ols_reg1 <- function(Y, X, L) {
  ok <- complete.cases(Y, X); Y <- Y[ok, , drop = FALSE]; X <- X[ok, , drop = FALSE]; N <- nrow(X); Kx <- ncol(X); qx <- qr(X, tol = 1e-13); iXX <- chol2inv(qr.R(qx))
  B <- qr.coef(qx, Y); err <- Y - X %*% B; se <- matrix(0, Kx, ncol(Y))
  for (v in 1:ncol(Y)) { errs <- X * err[, v]; V0 <- crossprod(errs) / N
    for (i in 1:L) { S <- crossprod(errs[1:(N - i), , drop = FALSE], errs[(1 + i):N, , drop = FALSE]) / N; V0 <- V0 + (1 - i / (L + 1)) * (S + t(S)) }
    se[, v] <- sqrt(diag(N * iXX %*% V0 %*% iXX)) }
  list(beta = B, se = se, err = err, N = N, K = Kx)
}
lagmat <- function(Z, ls) do.call(cbind, lapply(ls, function(l) rbind(matrix(NA, l, ncol(Z)), Z[1:(nrow(Z) - l), , drop = FALSE])))
lead <- function(Z, h) rbind(Z[(1 + h):nrow(Z), , drop = FALSE], matrix(NA, h, ncol(Z)))
# ---- L1: engine: varm lp (shock = gs1, contemp = logip, logcpi) vs ols_reg on the same design, all h, all responses ----
j <- 3; ct <- 1:2
Wd <- cbind(y[, j], 1, lagmat(y, 1:lags), y[, ct])                       # shock first; lp() orders [shock, const, lags, contemp]
res <- lp(y, shock = j, horizon = H, lags = lags, contemp = ct, se = "nw", nw_lag = NULL)  # varm default L = h+1
dB <- 0; dS <- 0; relmax <- function(a, b, se) { m <- se > 1e-5; if (any(m)) max(abs(a - b)[m] / se[m]) else 0 }
for (h in 0:H) { r <- ols_reg1(lead(y, h), Wd, h + 1); dB <- max(dB, relmax(r$beta[1, ], res$estimate[h + 1, ], res$se[h + 1, ])); dS <- max(dS, relmax(r$se[1, ], res$se[h + 1, ], res$se[h + 1, ])) }
rep_("L1a LP beta (shock gs1, contemp logip/logcpi, 4 lags) equals QR-OLS on the same design, h=0..12 (rel. to se, tol 1e-5)", dB < 1e-5, "max|diff|/se=", signif(dB, 3))
rep_("L1b LP Newey-West se (varm nw_lag = h+1) equals Hamilton se of ols_reg (rel., tol 1e-5)", dS < 1e-5, "max rel diff=", signif(dS, 3))
res2 <- lp(y, shock = j, horizon = H, lags = lags, contemp = ct, se = "nw", nw_lag = lags + 1)
dS2 <- 0; for (h in 0:H) { r <- ols_reg1(lead(y, h), Wd, lags + h + 1); res3 <- lp(y, shock = j, horizon = h, lags = lags, contemp = ct, se = "nw", nw_lag = lags + h + 1); dS2 <- max(dS2, relmax(r$se[1, ], res3$se[h + 1, ], res3$se[h + 1, ])) }
rep_("L1c same with BVAR_'s bandwidth L = lags + h + 1 (nw_lag argument; rel., tol 1e-5)", dS2 < 1e-5, "max rel diff=", signif(dS2, 3))
# ---- L2: BVAR_ 'LP-Cholesky' (ir_lp): beta_{h-1}(lag-1 block)' * chol(Sigma_0), vs varm unit-impact LP with contemporaneous controls ----
r0 <- ols_reg1(lead(y, 0), cbind(lagmat(y, 1:lags), 1), 1); S0 <- crossprod(r0$err) / (r0$N - r0$K); Om <- t(chol(S0))
irB <- array(0, c(K, H + 1, K)); irB[, 1, ] <- Om
for (h in 1:H) { rr <- ols_reg1(lead(y, h - 1), cbind(lagmat(y, 1:lags), 1), 1); irB[, h + 1, ] <- t(rr$beta[1:K, ]) %*% Om }
cmp <- function(j) { cj <- if (j > 1) 1:(j - 1) else NULL; l <- lp(y, shock = j, horizon = H, lags = lags, contemp = cj); v <- l$estimate   # (H+1) x K, unit impact of y_j
  bn <- t(irB[, , j]) / Om[j, j]; { a <- (abs(v[, j:K, drop = FALSE] - bn[, j:K, drop = FALSE]) / pmax(l$se[, j:K, drop = FALSE], 1e-8))[-1, , drop = FALSE]; max(a) } }
dd <- sapply(1:K, cmp)
fitv <- varm(y, lags); vi <- irf(fitv, id_chol(), horizon = H, ci = "none")$estimate
dv <- sapply(1:K, function(j) max(abs(t(vi[, j, ]) / vi[j, j, 1] - t(irB[, , j]) / Om[j, j])))
cat("   max|varm unit-impact LP - BVAR_ ir_lp/Omega_jj|/se per shock (h>=1):", signif(dd, 3), " ; BVAR_ ir_lp vs VAR-Cholesky IRF (varm) normalised:", signif(dv, 3), "\n")
rep_("L2a BVAR_ LP-Cholesky response vs varm LP (different regressions: lag-1 coefficient of y_{t+h-1} vs contemporaneous controls): differences < 1 se, h>=1", max(dd) < 1, "max|diff|/se_varm=", signif(max(dd), 3), " (sampling difference, not an algebraic identity)")
# ---- L3: proxy LP: BVAR_ regresses y_{t+h} on [lags, proxy, const]; varm lp(shock = numeric proxy) adds `lags` lags of the proxy ----
fd <- suppressMessages(as.data.frame(readxl::read_excel(file.path(EXDIR, "factor_data.xlsx"), sheet = "factor_data"))); prox <- fd$ff4_tc; Tm <- length(prox)
ys <- y[(T - Tm + 1):T, ]; yp <- ys[, c(3, 1, 2, 4)]    # policy first not required for LP
Xb <- function(h) { ok <- NULL; cbind(lagmat(ys, 1:lags), prox, 1) }
dBp <- 0; for (h in c(0, 3, 6, 12)) { rb <- ols_reg1(lead(ys, h), Xb(), lags + h + 1)
  vv <- lp(ys, shock = prox, horizon = h, lags = lags, se = "nw", nw_lag = lags + h + 1)
  dBp <- max(dBp, max(abs(rb$beta[K * lags + 1, ] - vv$estimate[h + 1, ]) / vv$se[h + 1, ])) }
cat("   proxy LP: max|beta_BVAR_ - beta_varm|/se over h in {0,3,6,12} and 4 responses:", signif(dBp, 3), " (varm also controls for", lags, "lags of the proxy)\n")
rep_("L3 proxy-LP coefficients agree within 0.5 se (difference = lags of the proxy included by varm only)", dBp < 0.5, "max|diff|/se=", signif(dBp, 3))
# engine identity with proxy lags: varm vs lm on the exact varm design
dE <- 0; for (h in c(0, 6, 12)) { Xv <- cbind(prox, 1, lagmat(ys, 1:lags), lagmat(matrix(prox), 1:lags)); rr <- ols_reg1(lead(ys, h), Xv, h + 1)
  vv <- lp(ys, shock = prox, horizon = h, lags = lags, se = "nw"); dE <- max(dE, max(abs(rr$beta[1, ] - vv$estimate[h + 1, ]) / vv$se[h + 1, ]), max(abs(rr$se[1, ] - vv$se[h + 1, ]) / vv$se[h + 1, ])) }
rep_("L3b varm proxy LP engine == QR-OLS/NW on its own design (beta and se, rel. 1e-5)", dE < 1e-5, "max rel diff=", signif(dE, 3))
# ---- L4: confidence band convention ----
cat("   bands: BVAR_ uses beta +- t_{0.95,T-ny-1} * se (talpha = tinv), varm uses normal quantiles: t=", signif(abs(qt(0.05, T - K - 1)), 5), " vs z=", signif(qnorm(.95), 5), "\n")
