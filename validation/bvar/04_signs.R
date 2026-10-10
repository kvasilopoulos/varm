# Sign restrictions: varm id_sign() vs BVAR_ iresponse_sign.m / generateQ.m / checkrestrictions2.m (example_3_irf.m case 2), at the OLS point estimate.
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm()
set.seed(20261011)
GK <- read_mat5(file.path(EXDIR, "DataGK.mat")); y <- cbind(logip = GK$logip, logcpi = GK$logcpi, gs1 = GK$gs1, ebp = GK$ebp); lags <- 12; K <- 4; H <- 11
o <- yxb(y, lags); Phi <- qr.coef(qr(o$X, tol = 1e-12), o$Y); u <- o$Y - o$X %*% Phi; Sig <- crossprod(u) / (nrow(o$Y) - ncol(o$X))
fit <- varm(y, lags); fit$B <- Phi; fit$U <- u; fit$Sigma <- Sig
L <- t(chol(Sig)); MA <- varm:::cpp_ma(Phi, K, lags, H)       # [K,K,H+1]
haar_B <- function() { q <- qr(matrix(rnorm(K * K), K)); Q <- qr.Q(q); Q %*% diag(sign(diag(qr.R(q)))) }   # generateQ.m
# BVAR_ check: restriction strings like y(var, hs, shock) > 0 / < 0; y is the K x hor0 x K IRF array. Global sign flip of the whole array (checkrestrictions2.m)
bvar_check <- function(Y, rs) {            # rs: list(var, hs (1-based, impact = 1), shock, sign)
  ok <- function(Y) sum(sapply(rs, function(r) { v <- Y[r$var, r$hs, r$shock] * r$sign; as.numeric(min(v > 0)) })) == length(rs)
  p <- ok(Y); m <- ok(-Y); c(accept = p || m, flip = if (m && !p) -1 else 1)
}
irfs_B <- function(Om, hor0) { A <- array(0, c(K, hor0, K)); for (h in 1:hor0) A[, h, ] <- MA[, , h] %*% L %*% Om; A }
run_bvar <- function(rs, N, hor0 = 3) {
  acc <- 0; keep <- list()
  for (i in 1:N) { Om <- haar_B(); r <- bvar_check(irfs_B(Om, hor0), rs); if (r["accept"] == 1) { acc <- acc + 1; keep[[acc]] <- r["flip"] * Om } }
  list(rate = acc / N, Om = keep)
}
mkS <- function(rs, hor0) { S <- array(0, c(K, K, hor0)); for (r in rs) S[r$var, r$shock, r$hs] <- r$sign; S }
# ---- S1: single-shock restrictions: gs1 up (h=1..3) and cpi down (h=1..3) after shock 1 (y(3,1:3,1)>0; y(2,1:3,1)<0) ------------------------------
rs1 <- list(list(var = 3, hs = 1:3, shock = 1, sign = 1), list(var = 2, hs = 1:3, shock = 1, sign = -1))
N <- 60000; rb <- run_bvar(rs1, N)
rv <- irf(fit, id_sign(mkS(rs1, 3), nrot = 6000, maxtry = 600000), horizon = H, ci = "none")
nacc <- dim(rv$draws)[4]
rate_v <- rv$accept_rate
se <- sqrt(rb$rate * (1 - rb$rate) / N)
cat(sprintf("   S1 acceptance: BVAR_ rule %.4f (+-%.4f)  varm %.4f ; accepted draws varm %d, BVAR_ %d\n", rb$rate, se, rate_v, nacc, length(rb$Om)))
rep_("S1a acceptance probability equal (single-shock restrictions)", abs(rb$rate - rate_v) < 4 * sqrt(se^2 + rate_v * (1 - rate_v) / 1e5), "diff=", signif(rate_v - rb$rate, 3))
# distribution of accepted IRFs of shock 1 (columns of Theta_h L Omega), all variables, horizons 0..H
IB <- sapply(rb$Om, function(Om) { a <- array(0, c(K, H + 1)); for (h in 0:H) a[, h + 1] <- (MA[, , h + 1] %*% L %*% Om)[, 1]; as.vector(a) })   # (K*(H+1)) x n
IV <- matrix(rv$draws[, 1, , ], K * (H + 1))
zz <- (rowMeans(IB) - rowMeans(IV)) / sqrt(apply(IB, 1, var) / ncol(IB) + apply(IV, 1, var) / ncol(IV))
qd <- sapply(c(.05, .16, .5, .84, .95), function(q) { a <- apply(IB, 1, quantile, q); b <- apply(IV, 1, quantile, q); max(abs(a - b) / pmax(apply(IB, 1, sd), 1e-12)) })
rep_("S1b accepted IRF (shock 1): means agree (|z|<4.5, 48 cells) and quantiles within 0.1 sd", max(abs(zz)) < 4.5 && max(qd) < 0.1, "max|z|=", round(max(abs(zz)), 2), " max quantile gap/sd=", signif(max(qd), 3))
# ---- S2: restrictions on TWO shocks: BVAR_ flips the sign of the whole array (global fsign), varm flips per column ---------------------------------
rs2 <- c(rs1, list(list(var = 4, hs = 1:2, shock = 2, sign = 1), list(var = 3, hs = 1:2, shock = 2, sign = -1)))
rb2 <- run_bvar(rs2, 60000)
rv2 <- irf(fit, id_sign(mkS(rs2, 3), nrot = 3000, maxtry = 600000), horizon = H, ci = "none")
cat(sprintf("   S2 acceptance (restrictions on 2 shocks): BVAR_ global-flip rule %.4f  varm per-column rule %.4f  ratio %.2f (theory 2 under Haar sign symmetry)\n", rb2$rate, rv2$accept_rate, rv2$accept_rate / rb2$rate))
rep_("S2 two-shock restrictions: acceptance rate ratio varm/BVAR_ ~ 2 (BVAR_ cannot flip columns independently)", abs(rv2$accept_rate / rb2$rate - 2) < 0.25, "ratio=", round(rv2$accept_rate / rb2$rate, 3))
# ---- S3: every accepted varm rotation is orthogonal & satisfies the restrictions -----------------------------------------------------------------
P <- rv$draw; Pm <- varm:::cpp_sign_irf(array(Phi, c(dim(Phi), 1)), array(Sig, c(K, K, 1)), lags, 2L, mkS(rs1, 3), 200L, 20000L)$P
okP <- all(sapply(1:dim(Pm)[3], function(i) mx(Pm[, , i] %*% t(Pm[, , i]), Sig) < 1e-10))
sg <- sapply(1:dim(Pm)[3], function(i) { A <- sapply(0:2, function(h) (MA[, , h + 1] %*% Pm[, , i])[, 1]); all(A[3, ] > 0) && all(A[2, ] < 0) })
rep_("S3 varm accepted impact matrices: P P' = Sigma and signs hold", okP && all(sg), "n=", dim(Pm)[3])
