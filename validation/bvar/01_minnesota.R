# Minnesota prior / marginal likelihood / posterior: varm::bvarm(prior_minnesota) vs BVAR_ (bvar_ml.m, varprior.m, rfvar3.m, bvar_.m)
# Run: Rscript validation/bvar/01_minnesota.R
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm()
set.seed(20261011)
D <- read_mat5(file.path(EXDIR, "Data.mat"))
y <- cbind(IPI = D$IPI, HICP = D$HICP, CORE = D$CORE)   # first three variables of the euro-area tutorial data (T = 300)
cat("data: Data.mat, T =", nrow(y), " K = 3\n")
# mapping BVAR_ -> varm:  lambda = sqrt(omega)/tau, alpha = 2*decay, soc = 1/mu, sur = 1/lambda_B, psi = omega*sig^2
# (sig = std of the first `lags` observations, as in bvar_ml.m).  See diary for the derivation.
vm <- function(y, lags, h, psi, w, var = 1e7)
  suppressWarnings(bvarm(y, p = lags, ndraw = 5, prior = prior_minnesota(lambda = sqrt(w) / h[1], alpha = 2 * h[2],
    soc = 1 / h[4], sur = 1 / h[3], psi = psi, hyper = character(0), var = var)))
w <- 2   # M = 3, omega = 2 is the only integer configuration whose inverse-Wishart dof coincide (M*omega = M + 3)
grid <- list(c(3, .5, 5, 2, w), c(1, .5, 5, 2, w), c(10, 1, 1, .5, w), c(.3, 2, 20, 5, w), c(5, 1, 1, .5, w))
# ---- T1: ML equals exact NIW ML of BVAR_'s own prior (sanity of the R port) ----------------------
lmg <- function(a, M) M * (M - 1) / 4 * log(pi) + sum(lgamma(a - (0:(M - 1)) / 2))
niw_ml <- function(Y, X, Bp, Op, Sp, nup) {
  N <- nrow(Y); M <- ncol(Y); Oi <- solve(Op); P <- crossprod(X) + Oi; Bh <- solve(P, crossprod(X, Y) + Oi %*% Bp)
  Sq <- Sp + crossprod(Y - X %*% Bh) + t(Bh - Bp) %*% Oi %*% (Bh - Bp); nuq <- nup + N
  ld <- function(A) as.numeric(determinant(A)$modulus)
  -M * N / 2 * log(pi) + M / 2 * (-ld(P) + ld(Oi)) + nup / 2 * ld(Sp) - nuq / 2 * ld(Sq) + lmg(nuq / 2, M) - lmg(nup / 2, M)
}
for (lags in c(2, 4)) {
  b <- bvar_ml(grid[[1]], y, lags); o <- yxb(y, lags)
  e <- niw_ml(o$Y, o$X, b$prior$B, b$prior$xxi, crossprod(b$prior$u), b$prior$df)
  rep_(sprintf("T1 BVAR_ ML == exact NIW ML of its dummy prior (p=%d)", lags), abs(e - b$logml) < 1e-6, "BVAR_=", format(b$logml, digits = 10), " exact=", format(e, digits = 10))
}
# ---- T2: varm logml vs BVAR_ logml over a hyper-parameter grid ----------------------------------
for (lags in c(2, 4)) for (h in grid) {
  b0 <- bvar_ml(h, y, lags); b1 <- bvar_ml(h, y, lags, ybar_shift = lags); f <- vm(y, lags, h, w * b0$sig^2, w)
  d0 <- f$logml - b0$logml; d1 <- f$logml - b1$logml
  cat(sprintf("   p=%d h(tau,decay,lambda,mu,omega)=%-18s BVAR_=%.4f varm=%.4f | diff(ybar=y[1:p])=%+.5f  diff(ybar=y[p+1:2p])=%+.2e\n",
              lags, paste(h, collapse = ","), b0$logml, f$logml, d0, d1))
  if (identical(h, grid[[1]]) && lags == 2) { T2a <- abs(d0); T2b <- abs(d1) }
}
rep_("T2a varm ML == BVAR_ ML (BVAR_ convention ybar = mean(y[1:p]))", T2a < 1e-4, "|diff| p=2 grid1 =", signif(T2a, 3), " (diary: nonzero, hyper-dependent, -0.007..-0.023)")
rep_("T2b varm ML == BVAR_ ML if ybar = mean(y[p+1:2p]) (what varm uses)", T2b < 1e-4, "|diff|=", signif(T2b, 3))
# ---- T3: dependence on the constant's prior variance --------------------------------------------
f7 <- vm(y, 2, grid[[1]], w * bvar_ml(grid[[1]], y, 2)$sig^2, w, 1e7)$logml; f12 <- vm(y, 2, grid[[1]], w * bvar_ml(grid[[1]], y, 2)$sig^2, w, 1e12)$logml
rep_("T3 varm logml insensitive to var (1e7 vs 1e12), i.e. matches BVAR_'s flat constant", abs(f7 - f12) < 1e-6, "diff=", signif(f7 - f12, 3))
# ---- T4: Sigma prior family (M = 2: dof of BVAR_'s omega-dummy prior cannot equal GLP's M+2 for the same scale) --
y2 <- y[, 1:2]; dd <- sapply(grid[1:4], function(h) { h[5] <- 1; b <- bvar_ml(h, y2, 2, ybar_shift = 2); f <- vm(y2, 2, h, b$sig^2, 1); f$logml - b$logml })
rep_("T4 (M=2, omega=1) ML differs hyper-dependently (different IW dof) -- convention difference", max(dd) - min(dd) > 1e-3, "diff range:", paste(signif(range(dd), 4), collapse = " .. "))
# ---- T5: posterior parameters (B-hat, IW scale, dof) at M=3, omega=2, ybar matched ---------------
lags <- 2; h <- grid[[1]]; b1 <- bvar_ml(h, y, lags, ybar_shift = lags); ols <- varm:::cpp_var_ols(y, lags, matrix(1, nrow(y), 1))
bb <- matrix(0, 7, 3); bb[cbind(1:3, 1:3)] <- 1
vv <- varm:::cpp_niw_fixed(c(sqrt(w) / h[1], 2 * h[2], 1 / h[4], 1 / h[3], w * b1$sig^2), ols$Y, ols$X, lags, 1L, TRUE, bb, 1e7, colMeans(ols$Y[seq_len(lags), , drop = FALSE]), 1, 5L)
rep_("T5a posterior mean B-hat equal", mx(vv$bhat, b1$post$B) < 1e-6, "max|diff|=", signif(mx(vv$bhat, b1$post$B), 3))
rep_("T5b IW scale S equal", mx(vv$S, crossprod(b1$post$u)) / max(abs(vv$S)) < 1e-6, "max rel diff=", signif(mx(vv$S, crossprod(b1$post$u)) / max(abs(vv$S)), 3))
rep_("T5c posterior IW dof equal (varm nu vs BVAR_ posterior.df)", vv$nu == b1$post$df, "varm nu=", vv$nu, " BVAR_ df=", b1$post$df)
# ---- T6: posterior draws: varm sampler vs BVAR_ sampler (bvar_.m draw loop) -----------------------
N <- 20000
bvar_draw <- function(post, N) {          # rand_inverse_wishart + Phi = PhiHat + chol(XXi)' Z chol(Sigma)'
  M <- ncol(post$B); k <- nrow(post$B); Sinv_chol <- chol(solve(crossprod(post$u))); XL <- t(chol(post$xxi))
  Bs <- array(0, c(k, M, N)); Ss <- array(0, c(M, M, N))
  for (n in 1:N) {
    X <- matrix(rnorm(post$df * M), post$df, M) %*% Sinv_chol; s <- svd(X); G <- (s$v * rep(1 / s$d^2, each = M)) %*% t(s$v)
    Ss[, , n] <- G; Bs[, , n] <- post$B + XL %*% matrix(rnorm(k * M), k, M) %*% chol(G)
  }
  list(B = Bs, Sigma = Ss)
}
fv <- suppressWarnings(bvarm(y, p = lags, ndraw = N, prior = prior_minnesota(lambda = sqrt(w) / h[1], alpha = 2 * h[2], soc = 1 / h[4], sur = 1 / h[3],
        psi = w * b1$sig^2, hyper = character(0))))
mb <- bvar_draw(b1$post, N)
flat <- function(a) t(matrix(a, prod(dim(a)[1:2])))
zf <- function(a, b) (colMeans(a) - colMeans(b)) / sqrt(apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b))
zB <- zf(flat(fv$draws$B), flat(mb$B)); zS <- zf(flat(fv$draws$Sigma), flat(mb$Sigma))
vrB <- apply(flat(fv$draws$B), 2, var) / apply(flat(mb$B), 2, var)
rep_("T6 draw moments (means of B, Sigma) agree: |z| < 4.5", max(abs(c(zB, zS))) < 4.5, "max|z| B=", round(max(abs(zB)), 2), " Sigma=", round(max(abs(zS)), 2), " var-ratio range", paste(round(range(vrB), 3), collapse = ".."))
# ---- T7: ML maximisation over (tau, lambda) with decay = 1, mu = 2 fixed: argmax of the two ML surfaces ----
sg <- bvar_ml(c(1, 1, 1, 1, w), y, 2)$sig
objB <- function(x, shift) -bvar_ml(c(exp(x[1]), 1, exp(x[2]), 2, w), y, 2, ybar_shift = shift)$logml
objV <- function(x) -vm(y, 2, c(exp(x[1]), 1, exp(x[2]), 2, w), w * sg^2, w)$logml
x0 <- log(c(3, 5)); ctl <- list(reltol = 1e-12, maxit = 2000)
oB1 <- optim(x0, objB, shift = 2, method = "Nelder-Mead", control = ctl); oB0 <- optim(x0, objB, shift = 0, method = "Nelder-Mead", control = ctl)
oV <- optim(x0, objV, method = "Nelder-Mead", control = ctl)
cat("   argmax (tau, lambda): BVAR_(ybar=y[p+1:2p])=", round(exp(oB1$par), 3), " varm=", round(exp(oV$par), 3), " BVAR_(ybar=y[1:p])=", round(exp(oB0$par), 3), "\n")
cat("   max logML: BVAR_ shifted", -oB1$value, " varm", -oV$value, " BVAR_ own", -oB0$value, "\n")
rep_("T7a ML maximisers agree when ybar matched (rel 1%)", max(abs(exp(oB1$par) / exp(oV$par) - 1)) < 0.01, "")
rep_("T7b ML maximiser with BVAR_'s own ybar convention differs from varm's by < 1%", max(abs(exp(oB0$par) / exp(oV$par) - 1)) < 0.01, "max rel diff", signif(max(abs(exp(oB0$par) / exp(oV$par) - 1)), 3))
# ---- T8: trending data (DataGLP: log real per-capita GDP, consumption, investment; p = 5, BVAR_ defaults tau=3, decay=.5, lambda=5, mu=2, omega=2) --
G <- read_mat5(file.path(EXDIR, "DataGLP.mat")); yg <- cbind(Y = G$Y, C = G$C, I = G$I); h <- c(3, .5, 5, 2, 2)
b0 <- bvar_ml(h, yg, 5); b1 <- bvar_ml(h, yg, 5, ybar_shift = 5); f <- vm(yg, 5, h, 2 * b0$sig^2, 2)
cat(sprintf("   BVAR_=%.4f varm=%.4f diff(BVAR_ ybar)=%+.4f  diff(varm ybar)=%+.2e ; ybar y[1:5]=%s  y[6:10]=%s\n", b0$logml, f$logml, f$logml - b0$logml, f$logml - b1$logml,
            paste(round(colMeans(yg[1:5, ]), 3), collapse = ","), paste(round(colMeans(yg[6:10, ]), 3), collapse = ",")))
rep_("T8a varm ML == BVAR_ ML with BVAR_'s ybar convention", abs(f$logml - b0$logml) < 1e-4, "diff=", signif(f$logml - b0$logml, 3))
rep_("T8b varm ML == BVAR_ ML with varm's ybar rows (confirms cause)", abs(f$logml - b1$logml) < 1e-4, "diff=", signif(f$logml - b1$logml, 3))
# ---- T9: sum-of-coefficients switched off (BVAR_ mu = 0  <->  varm soc = NULL), M=3, omega=2 (dof coincide: M*omega = M+3 again) ----
# Without the soc dummies nothing but the sur row pins down the constant direction, so the constant's prior variance `var` matters:
# BVAR_ treats the constant as flat; varm's default var = 1e7 is not flat enough here.
t9 <- function(v) sapply(grid[1:4], function(h) { h[4] <- 0; b <- bvar_ml(h, y, 2, ybar_shift = 2)
  f <- suppressWarnings(bvarm(y, p = 2, ndraw = 5, prior = prior_minnesota(lambda = sqrt(w) / h[1], alpha = 2 * h[2], soc = NULL, sur = 1 / h[3], psi = w * b$sig^2, hyper = character(0), var = v)))
  f$logml - b$logml })
d9a <- t9(1e7); d9b <- t9(1e15)
rep_("T9a soc = NULL, varm default var = 1e7 vs BVAR_ mu = 0: ML agree (1e-3)", max(abs(d9a)) < 1e-3, "diffs:", paste(signif(d9a, 3), collapse = ", "), " (varm ML depends on var when soc is off)")
rep_("T9b soc = NULL with var = 1e15 (flat constant) vs BVAR_ mu = 0 (ybar matched)", max(abs(d9b)) < 1e-4, "diffs:", paste(signif(d9b, 3), collapse = ", "))
# ---- T10: posterior Cholesky IRFs (h = 0, 4, 12, 24): varm draws vs BVAR_ sampler draws (M=3, omega=2, ybar matched), posterior means ----
ii <- sample(N, 4000)
irB <- sapply(ii, function(n) { a <- iresponse(mb$B[1:6, , n], mb$Sigma[, , n], 25, diag(3)); as.vector(a[, c(1, 5, 13, 25), ]) })
Bv <- fv$draws$B; Sv <- fv$draws$Sigma
irV <- sapply(ii, function(n) { a <- iresponse(Bv[1:6, , n], Sv[, , n], 25, diag(3)); as.vector(a[, c(1, 5, 13, 25), ]) })
zi <- (rowMeans(irB) - rowMeans(irV)) / sqrt(apply(irB, 1, var) / length(ii) + apply(irV, 1, var) / length(ii))
rep_("T10 posterior mean of Cholesky IRFs (36 cells) agree: max|z| < 4.5", max(abs(zi[is.finite(zi)])) < 4.5, "max|z|=", round(max(abs(zi[is.finite(zi)])), 2))
