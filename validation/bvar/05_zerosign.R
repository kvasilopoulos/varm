# Zero + sign restrictions (impact): varm id_zerosign() vs BVAR_ iresponse_zeros_signs.m (Binning's algorithm: findQs/findP/generateDraw),
# restriction set of example_3_irf.m case 4 on Gertler-Karadi data, at the OLS point estimate.
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm()
set.seed(20261011)
GK <- read_mat5(file.path(EXDIR, "DataGK.mat")); y <- cbind(logip = GK$logip, logcpi = GK$logcpi, gs1 = GK$gs1, ebp = GK$ebp); lags <- 12; K <- 4
o <- yxb(y, lags); Phi <- qr.coef(qr(o$X, tol = 1e-12), o$Y); u <- o$Y - o$X %*% Phi; Sig <- crossprod(u) / (nrow(o$Y) - ncol(o$X))
C1 <- t(chol(Sig)); beta <- Reduce(`+`, lapply(1:lags, function(l) t(Phi[(l - 1) * K + 1:K, ])))
# sign matrix as in sign2matrix: sr (K x K, NaN = free), ys (1 = no short-run zero), yr (1 = no long-run zero)
sr <- matrix(NA, K, K); sr[1, 1] <- 1; sr[2, 1] <- 1; sr[3, 1] <- 1; sr[1, 2] <- 1; sr[2, 2] <- -1; sr[3, 3] <- 1; sr[4, 3] <- 1
ys <- matrix(1, K, K); ys[1, 3] <- 0; ys[2, 3] <- 0; yr <- matrix(1, K, K); f <- rbind(ys, yr)
# ---- findQs.m ----
findQs <- function(k, f) {
  Qi <- lapply(1:k, function(i) { m <- diag(as.numeric(f[, i] == 0), nrow(f)); m[rowSums(m) > 0, , drop = FALSE] })
  rk <- sapply(Qi, function(m) if (nrow(m)) qr(m)$rank else 0)
  ord <- order(rk, decreasing = TRUE, method = "radix"); new <- rk[ord]
  flag <- if (any(new - (k - (1:k)) > 0)) 1 else if (all(new - (k - (1:k)) == 0)) 0 else -1
  index <- integer(k); for (i in 1:k) index[ord[i]] <- i
  list(Q = Qi[ord], index = index, flag = flag)
}
# ---- findP.m / generateDraw.m ----
nullcol <- function(M, k) { if (is.null(M) || nrow(M) == 0) return(c(1, rep(0, k - 1))); qq <- qr.Q(qr(t(M)), complete = TRUE); qq[, k] }
findP <- function(C, Qs, k) {
  Linf <- solve(diag(k) - beta, C); F <- rbind(C, Linf); P <- matrix(0, k, k)
  for (ii in 1:k) { Qt <- if (ii == 1) Qs$Q[[ii]] %*% F else rbind(Qs$Q[[ii]] %*% F, t(P[, seq_len(ii - 1), drop = FALSE])); P[, ii] <- nullcol(Qt, k) }
  P[, Qs$index]
}
genDraw <- function(C1) { q <- qr(matrix(rnorm(K * K), K)); Q <- qr.Q(q); Q <- Q %*% diag(ifelse(diag(qr.R(q)) < 0, -1, 1)); C1 %*% Q }
Qs <- findQs(K, f); cat("   findQs flag (0 = exactly identified, -1 = under-identified, 1 = over):", Qs$flag, " zero-restriction rows per shock:", sapply(Qs$Q, nrow), " index:", Qs$index, "\n")
bvar_draw <- function() {
  C <- genDraw(C1); P <- findP(C, Qs, K); W <- C %*% P
  ok <- all(sapply(1:K, function(j) { s <- !is.na(sr[, j]); all(sign(W[s, j]) == sr[s, j]) }))
  list(ok = ok, W = W)
}
N <- 40000; res <- replicate(N, bvar_draw(), simplify = FALSE); okB <- sapply(res, `[[`, "ok"); WB <- res[okB]
cat(sprintf("   BVAR_ acceptance %.4f (%d of %d); note the sign of W columns is NOT normalised by BVAR_ (no flip), only accepted if signs match\n", mean(okB), sum(okB), N))
zB <- max(sapply(WB, function(r) max(abs(r$W[1:2, 3])))); pB <- max(sapply(WB, function(r) mx(r$W %*% t(r$W), Sig)))
rep_("Z1 BVAR_ draws: zeros exact and W W' = Sigma", zB < 1e-12 && pB < 1e-10, "max|zero|=", signif(zB, 3), " max|WW'-Sigma|=", signif(pB, 3))
# ---- varm ----
Sg <- array(0, c(K, K, 1)); Sg[1, 1, 1] <- 1; Sg[2, 1, 1] <- 1; Sg[3, 1, 1] <- 1; Sg[1, 2, 1] <- 1; Sg[2, 2, 1] <- -1; Sg[3, 3, 1] <- 1; Sg[4, 3, 1] <- 1
Zm <- matrix(0, K, K); Zm[1, 3] <- 1; Zm[2, 3] <- 1
id <- id_zerosign(sign = Sg, zero = Zm, nrot = 20000, maxtry = 2e6)
sv <- varm:::structural_irf(array(Phi, c(dim(Phi), 1)), array(Sig, c(K, K, 1)), K, lags, 0L, id, colnames(y))
WV <- sv$P; nV <- dim(WV)[3]
zV <- max(sapply(1:nV, function(i) max(abs(WV[1:2, 3, i])))); pV <- max(sapply(1:nV, function(i) mx(WV[, , i] %*% t(WV[, , i]), Sig)))
rep_("Z2 varm draws: zeros exact and P P' = Sigma", zV < 1e-12 && pV < 1e-10, "max|zero|=", signif(zV, 3), " max|PP'-Sigma|=", signif(pV, 3), " n=", nV)
cat(sprintf("   varm acceptance %.4f (per-column sign flip allowed => BVAR_ x 2^(#restricted shocks) expected: %.4f)\n", sv$accept_rate, mean(okB) * 2^3))
# ---- distributions of the 16 impact elements ----
MB <- sapply(WB, function(r) as.vector(r$W)); MV <- matrix(WV, K * K)
zz <- (rowMeans(MB) - rowMeans(MV)) / sqrt(apply(MB, 1, var) / ncol(MB) + apply(MV, 1, var) / ncol(MV))
qd <- sapply(c(.05, .25, .5, .75, .95), function(q) max(abs(apply(MB, 1, quantile, q) - apply(MV, 1, quantile, q)) / pmax(apply(MB, 1, sd), 1e-9)))
print(round(matrix(zz, K), 2))
rep_("Z3 distribution of accepted impact matrices (16 elements): |z| < 4.5 and quantile gaps < 0.1 sd", max(abs(zz[is.finite(zz)])) < 4.5 && max(qd) < 0.1, "max|z|=", round(max(abs(zz[is.finite(zz)])), 2), " max quantile gap/sd=", signif(max(qd), 3))
# Z3b: BVAR_ draws with per-column sign flipping before the check (isolates the sampling algorithm from sign handling)
bvar_draw_flip <- function() {
  C <- genDraw(C1); P <- findP(C, Qs0, K); W <- C %*% P
  for (j in 1:K) { s <- !is.na(sr0[, j]); if (any(s) && all(sign(W[s, j]) == -sr0[s, j])) W[, j] <- -W[, j] }
  ok <- all(sapply(1:K, function(j) { s <- !is.na(sr0[, j]); all(sign(W[s, j]) == sr0[s, j]) })); list(ok = ok, W = W)
}
Qs0 <- Qs; sr0 <- sr
resF <- replicate(N, bvar_draw_flip(), simplify = FALSE); okF <- sapply(resF, `[[`, "ok"); WF <- resF[okF]; MF <- sapply(WF, function(r) as.vector(r$W))
zzF <- (rowMeans(MF) - rowMeans(MV)) / sqrt(apply(MF, 1, var) / ncol(MF) + apply(MV, 1, var) / ncol(MV))
cat(sprintf("   BVAR_ with column flipping: acceptance %.4f (varm %.4f)\n", mean(okF), sv$accept_rate)); print(round(matrix(zzF, K), 2))
qdF <- sapply(c(.05, .25, .5, .75, .95), function(q) max(abs(apply(MF, 1, quantile, q) - apply(MV, 1, quantile, q)) / pmax(apply(MF, 1, sd), 1e-9)))
rep_("Z3b BVAR_ sampler (with column flips) vs varm: acceptance equal", abs(mean(okF) - sv$accept_rate) < 4 * sqrt(mean(okF) * (1 - mean(okF)) / N + 0.1 * 0.9 / 20000), "BVAR_flip=", signif(mean(okF), 4), " varm=", signif(sv$accept_rate, 4))
rep_("Z3c same: distribution of impact matrices |z|<4.5 and quantile gap < 0.1 sd", max(abs(zzF[is.finite(zzF)])) < 4.5 && max(qdF) < 0.1, "max|z|=", round(max(abs(zzF[is.finite(zzF)])), 2), " max q-gap/sd=", signif(max(qdF), 3))
# ---- Z4: exactly identified zero patterns (unique up to column sign): both implementations must reproduce the closed-form answer ----
colsign <- function(A, B) { s <- sign(colSums(A * B)); A %*% diag(s) }       # align column signs
fU <- rbind(1 - upper.tri(diag(K)), matrix(1, K, K))                          # impact zeros above the diagonal -> W lower triangular
sr0 <- matrix(NA, K, K); old <- list(sr = sr, f = f, Qs = Qs); sr <- sr0; f <- fU; Qs <- findQs(K, f)
Wb <- { C <- genDraw(C1); findP(C, Qs, K) -> P; C %*% P }
Wb <- colsign(Wb, C1)
rep_("Z4a BVAR_ (exactly identified impact zeros) reproduces Cholesky(Sigma) up to column sign", mx(Wb, C1) < 1e-10, "max|diff|=", signif(mx(Wb, C1), 3))
Zi <- matrix(NA, K, K); Zi[upper.tri(Zi)] <- 0
idr <- id_restrictions(impact = Zi)
Pv <- varm:::structural_irf(array(Phi, c(dim(Phi), 1)), array(Sig, c(K, K, 1)), K, lags, 0L, idr, colnames(y))$irf[, , 1, 1]
rep_("Z4b varm id_restrictions (impact zeros) reproduces Cholesky(Sigma)", mx(colsign(Pv, C1), C1) < 1e-10, "max|diff|=", signif(mx(colsign(Pv, C1), C1), 3))
# long-run zeros above the diagonal: Blanchard-Quah.  BVAR_ via findQs/findP vs closed form F*chol(C1 S C1')
fL <- rbind(matrix(1, K, K), 1 - upper.tri(diag(K))); f <- fL; Qs <- findQs(K, f)
Wl <- colsign({ C <- genDraw(C1); findP(C, Qs, K) -> P; C %*% P }, diag(K) - beta %*% diag(K))
Cinf <- solve(diag(K) - beta); BQ <- (diag(K) - beta) %*% t(chol(Cinf %*% Sig %*% t(Cinf)))
rep_("Z5a BVAR_ (exactly identified long-run zeros) reproduces Blanchard-Quah up to column sign", mx(colsign(Wl, BQ), BQ) < 1e-8, "max|diff|=", signif(mx(colsign(Wl, BQ), BQ), 3))
Li <- matrix(NA, K, K); Li[upper.tri(Li)] <- 0; idl <- id_restrictions(longrun = Li)
Pl <- varm:::structural_irf(array(Phi, c(dim(Phi), 1)), array(Sig, c(K, K, 1)), K, lags, 0L, idl, colnames(y))$irf[, , 1, 1]
rep_("Z5b varm id_restrictions (long-run zeros) reproduces Blanchard-Quah", mx(colsign(Pl, BQ), BQ) < 1e-8, "max|diff|=", signif(mx(colsign(Pl, BQ), BQ), 3))
# ---- Z6: unfiltered law of the zero-restricted column: q3 = L^{-1} W[,3] lives on a circle (2-dim plane); is its angle uniform? ----------
Qs <- Qs0; f <- rbind(ys, yr)
Bpl <- qr.Q(qr(solve(C1, diag(K)[, 3:4])))                      # orthonormal basis of the plane L^-1 span(e3,e4)
ang <- function(W3) { q <- solve(C1, W3); atan2(sum(q * Bpl[, 2]), sum(q * Bpl[, 1])) }
aB <- replicate(20000, { C <- genDraw(C1); P <- findP(C, Qs, K); ang((C %*% P)[, 3]) })
Wv <- varm:::cpp_zs_irf(array(Phi, c(dim(Phi), 1)), array(Sig, c(K, K, 1)), lags, 0L, array(0, c(K, K, 1)), array(Zm, c(K, K, 1)), matrix(0, K, K), 20000L, 20000L)$P
aV <- sapply(1:dim(Wv)[3], function(i) ang(Wv[, 3, i]))
fold <- function(a) (a %% pi) / pi                              # sign of q3 irrelevant -> fold to [0,1)
kB <- suppressWarnings(ks.test(fold(aB), "punif")); kV <- suppressWarnings(ks.test(fold(aV), "punif"))
rep_("Z6a varm: angle of the zero-restricted q3 within its null-space plane is uniform (KS)", kV$p.value > 0.001, "KS p=", signif(kV$p.value, 3), " D=", signif(kV$statistic, 3))
rep_("Z6b BVAR_ (findP with random C): same angle is uniform (KS)", kB$p.value > 0.001, "KS p=", signif(kB$p.value, 3), " D=", signif(kB$statistic, 3), "(non-uniform => different implied prior over rotations; neither sampler carries the ARRW importance weights)")
# ---- Z7: law of the NEXT column given q3 (conditional Haar => q1 uniform on the 2-sphere in q3's orthogonal complement => coordinate in a fixed basis is U(-1,1)) ----
coord <- function(Q) { q3 <- Q[, 3]; Bc <- qr.Q(qr(q3), complete = TRUE)[, 2:4]; sum(Bc[, 1] * Q[, 1]) }
QB <- replicate(20000, { C <- genDraw(C1); P <- findP(C, Qs, K); W <- C %*% P; coord(solve(C1, W)) })
QV <- sapply(1:dim(Wv)[3], function(i) coord(solve(C1, Wv[, , i])))
kB <- suppressWarnings(ks.test((QB + 1) / 2, "punif")); kV <- suppressWarnings(ks.test((QV + 1) / 2, "punif"))
rep_("Z7a varm: q1 | q3 has the conditional-Haar marginal (coordinate ~ U(-1,1)) (KS)", kV$p.value > 0.001, "KS p=", signif(kV$p.value, 3), " D=", signif(kV$statistic, 3))
rep_("Z7b BVAR_: q1 | q3 has the conditional-Haar marginal (KS)", kB$p.value > 0.001, "KS p=", signif(kB$p.value, 3), " D=", signif(kB$statistic, 3))
