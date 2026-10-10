# Faithful R ports of the deterministic parts of Cesa-Bianchi's VAR-Toolbox (VAR/*.m).
# Arrays follow MATLAB conventions: IR[h, var, shock] (h=1 is impact), VD in percent.

ml_makexy <- function(D, p, const) {            # VARmakexy.m: X = [const, trend, lag1..lagp]
  n <- nrow(D); Y <- D[(p + 1):n, , drop = FALSE]
  X <- do.call(cbind, lapply(1:p, function(l) D[(p + 1 - l):(n - l), , drop = FALSE]))
  if (const == 1) X <- cbind(1, X)
  if (const == 2) X <- cbind(1, seq_len(n - p), X)   # trend = 1..Te (NOT 1..T)
  list(Y = Y, X = X)
}
ml_var <- function(D, p, const = 1) {          # VARmodel.m (reduced form)
  xy <- ml_makexy(D, p, const); X <- xy$X; Y <- xy$Y
  Te <- nrow(Y); K <- ncol(D); ntot <- ncol(X)
  Ft <- solve(crossprod(X), crossprod(X, Y)); res <- Y - X %*% Ft
  F <- t(Ft); sigma <- crossprod(res) / (Te - ntot)
  Fcomp <- rbind(F[, (1 + const):(K * p + const), drop = FALSE],
                 if (p > 1) cbind(diag(K * (p - 1)), matrix(0, K * (p - 1), K)))
  list(D = D, p = p, K = K, const = const, X = X, Y = Y, Ft = Ft, F = F, sigma = sigma, resid = res,
       Fcomp = Fcomp, maxEig = max(Mod(eigen(Fcomp, only.values = TRUE)$values)), nobs = Te, ntot = ntot)
}
ml_wold <- function(v, nsteps) {               # compute_wold.m
  K <- v$K; p <- v$p; Fp <- array(0, c(K, K, nsteps))
  for (i in seq_len(min(p, nsteps))) Fp[, , i] <- v$F[, v$const + (i - 1) * K + 1:K]
  PSI <- array(0, c(K, K, nsteps)); PSI[, , 1] <- diag(K)
  for (i in seq_len(nsteps)[-1]) { a <- 0; for (j in 1:(i - 1)) a <- a + PSI[, , i - j] %*% Fp[, , j]; PSI[, , i] <- a }
  PSI
}
ml_ir <- function(v, B, nsteps, impact = 0, recurs = "wold") {   # compute_IR.m
  K <- v$K; PSI <- ml_wold(v, nsteps); IR <- array(NA, c(nsteps, K, K))
  for (m in 1:K) {
    imp <- numeric(K); imp[m] <- if (impact == 0) 1 else 1 / B[m, m]
    for (k in 1:nsteps) {
      if (recurs == "wold") IR[k, , m] <- PSI[, , k] %*% B %*% imp
      else { Fn <- diag(nrow(v$Fcomp)); for (i in seq_len(k - 1)) Fn <- Fn %*% v$Fcomp; IR[k, , m] <- (Fn[1:K, 1:K] %*% B %*% imp) }
    }
  }
  IR
}
ml_vd <- function(v, B, nsteps) {              # compute_VD.m (percent)
  K <- v$K; PSI <- ml_wold(v, nsteps); S <- v$sigma
  MSE <- array(0, c(K, K, nsteps)); MSE[, , 1] <- S
  for (n in seq_len(nsteps)[-1]) MSE[, , n] <- MSE[, , n - 1] + PSI[, , n] %*% S %*% t(PSI[, , n])
  VD <- array(0, c(nsteps, K, K))
  for (i in 1:K) {
    Ms <- B[, i] %o% B[, i]; acc <- Ms
    for (n in 1:nsteps) {
      if (n > 1) acc <- acc + PSI[, , n] %*% Ms %*% t(PSI[, , n])
      VD[n, i, ] <- 100 * diag(acc) / diag(MSE[, , n])
    }
  }
  VD
}
ml_B <- function(v, ident) {                   # recover_B.m short/long
  if (ident == "short") return(t(chol(v$sigma)))
  K <- v$K; Finf <- solve(diag(nrow(v$Fcomp)) - v$Fcomp)[1:K, 1:K]
  D <- t(chol(Finf %*% v$sigma %*% t(Finf))); solve(Finf, D)
}
ml_hd <- function(v, B) {                      # compute_HD.m: returns HD$shock[time, shock, var] (NaN-padded) and init/const
  K <- v$K; p <- v$p; nobs <- v$nobs; const <- v$const; Fc <- v$Fcomp
  eps <- solve(B, t(v$resid)); Icomp <- cbind(diag(K), matrix(0, K, (p - 1) * K))
  Bbig <- rbind(B, matrix(0, K * (p - 1), K))
  shock <- array(0, c(nobs + p, K, K))
  for (j in 1:K) {
    st <- matrix(0, K * p, nobs + 1)
    for (i in 2:(nobs + 1)) {
      e <- numeric(K); e[j] <- eps[j, i - 1]
      st[, i] <- Bbig %*% e + Fc %*% st[, i - 1]
      shock[p + i - 1, j, ] <- Icomp %*% st[, i]
    }
  }
  shock[seq_len(p), , ] <- NA
  Xl <- v$X[, (1 + const):(K * p + const)]
  init <- matrix(0, K * p, nobs + 1); init[, 1] <- Xl[1, ]
  for (i in 2:(nobs + 1)) init[, i] <- Fc %*% init[, i - 1]
  cc <- matrix(0, K * p, nobs + 1); CC <- c(v$F[, 1], numeric(K * (p - 1)))
  for (i in 2:(nobs + 1)) cc[, i] <- CC + Fc %*% cc[, i - 1]
  list(shock = shock, init = t(Icomp %*% init)[-1, ], const = t(Icomp %*% cc)[-1, ])
}
# ---- sign restrictions: SignRestrictions.m, ds shocks, one fixed (v) reduced form
ml_orthnorm <- function(n) { q <- qr(matrix(rnorm(n * n), n)); Q <- qr.Q(q); R <- qr.R(q); Q %*% diag(ifelse(diag(R) < 0, -1, 1), n) }
ml_sign_once <- function(v, SIGN, sr_hor, sr_rot = 500) {
  dy <- nrow(SIGN); ds <- ncol(SIGN); C <- t(chol(v$sigma)); tried <- 0
  Phi <- if (sr_hor > 1) ml_wold(v, sr_hor) else NULL
  repeat {
    if (tried >= sr_rot) return(list(B = NULL, tried = tried))
    tried <- tried + 1
    terma <- C %*% ml_orthnorm(dy); termaa <- terma; order <- 1:dy; TermA <- 0
    for (ii in 1:ds) for (jj in 1:dy) if (is.finite(terma[1, jj])) {
      CH <- if (sr_hor > 1) sapply(1:sr_hor, function(h) Phi[, , h] %*% termaa[, jj] * 1) else terma[, jj]
      CH <- as.matrix(CH)
      # MATLAB uses sign check on aux_irf (VAR.B = termaa at the time); replicate: sign of IRF of column jj
      if (sum(CH * SIGN[, ii] < 0) == 0) { TermA <- TermA + 1; order[ii] <- jj; terma[, jj] <- NA; break }
      else if (sum(-CH * SIGN[, ii] < 0) == 0) { TermA <- TermA + 1; terma[, jj] <- NA; termaa[, jj] <- -termaa[, jj]; order[ii] <- jj; break }
    }
    if (TermA == ds) return(list(B = termaa[, order], tried = tried))
  }
}
# ---- proxy / IV: recover_B.m 'iv' (z may contain NA; rows dropped like MATLAB)
ml_iv <- function(v, IV) {
  p <- v$p; K <- v$K; up <- v$resid[, 1]; uq <- v$resid[, -1, drop = FALSE]
  Zt <- as.matrix(IV)[-seq_len(p), , drop = FALSE]
  keep <- !is.na(Zt[, 1]); pp <- up[keep]; q <- uq[keep, , drop = FALSE]; Z <- Zt[keep, 1]
  fs <- lm(pp ~ Z); p_hat <- fitted(fs)
  sqsp <- sapply(seq_len(K - 1), function(i) coef(lm(q[, i] ~ p_hat))[2])
  Biv <- c(1, sqsp)
  pq <- cbind(pp, q); cen <- sweep(pq, 2, colMeans(pq)); sb <- crossprod(cen) / (nrow(pq) - v$ntot)
  S11 <- sb[1, 1]; S21 <- sb[-1, 1]; S22 <- sb[-1, -1]
  Q <- sqsp %o% sqsp * S11 - (S21 %o% sqsp + sqsp %o% S21) + S22
  d <- S21 - sqsp * S11
  sp <- sqrt(S11 - drop(t(d) %*% solve(Q, d)))
  Biv * sp * sign(coef(fs)[2])
}
# ---- bootstrap, VARmodel.m 'bs' / 'wild' (no centring, no dof rescale), returns IR draws [h,var,shock,draw]
ml_boot <- function(v, ident, nsteps, ndraws, method = "bs", impact = 0, scale = FALSE) {
  K <- v$K; p <- v$p; Te <- v$nobs; const <- v$const; out <- array(NA, c(nsteps, K, K, ndraws)); tt <- 1
  res0 <- if (scale) v$resid * sqrt(Te / (Te - v$ntot)) else v$resid   # scale=TRUE mimics varm (not MATLAB)
  while (tt <= ndraws) {
    u <- if (method == "bs") res0[sample.int(Te, Te, TRUE), , drop = FALSE] else res0 * sample(c(-1, 1), Te, TRUE)
    ya <- matrix(0, Te + p, K); ya[1:p, ] <- v$D[1:p, ]
    for (j in (p + 1):(Te + p)) {
      lag <- as.vector(t(ya[(j - 1):(j - p), , drop = FALSE])); x <- c(if (const > 0) 1, if (const > 1) j - p, lag)
      ya[j, ] <- x %*% v$Ft + u[j - p, ]
    }
    vd <- ml_var(ya, p, const)
    if (vd$maxEig < .9999) { out[, , , tt] <- ml_ir(vd, ml_B(vd, ident), nsteps, impact); tt <- tt + 1 }
  }
  out
}
ml_prctile <- function(x, pr) quantile(x, pr, type = 5, names = FALSE)   # MATLAB prctile = Hyndman-Fan type 5
# ---- LP, LPmodel.m lp_single (OLS and IV, no longdiff), nlag lags of CTRL, impact=1
ml_lp <- function(endo, treat, ctrl, nlag, H, iv = NULL, const = 1) {
  n <- length(endo); endo_t <- endo[-seq_len(nlag)]; s <- treat[-seq_len(nlag)]
  cl <- do.call(cbind, lapply(1:nlag, function(l) ctrl[(nlag + 1 - l):(n - l), , drop = FALSE]))
  C <- cbind(if (const >= 1) 1, cl); nT <- length(endo_t)
  nw <- function(X, e, L) { hh <- X * e; G <- crossprod(hh); for (a in seq_len(L)) { za <- crossprod(hh[(a + 1):nrow(hh), , drop = FALSE], hh[1:(nrow(hh) - a), , drop = FALSE]); G <- G + (L + 1 - a) / (L + 1) * (za + t(za)) }; xi <- solve(crossprod(X)); sqrt(diag(xi %*% G %*% xi)) }
  res <- matrix(NA, H, 2)
  for (h in 1:H) {
    Y <- endo_t[h:nT]; nh <- length(Y)
    if (is.null(iv)) { X <- cbind(s[1:nh], C[1:nh, , drop = FALSE]); b <- solve(crossprod(X), crossprod(X, Y)); e <- Y - X %*% b
      res[h, ] <- c(b[1], nw(X, drop(e), h - 1)[1]) }
    else { Ch <- C[1:nh, , drop = FALSE]; rs <- function(a) a - Ch %*% solve(crossprod(Ch), crossprod(Ch, a))
      ry <- rs(Y); rd <- rs(s[1:nh]); z <- iv[-seq_len(nlag)][1:nh]; rz <- rs(z)
      dh <- rz * drop(crossprod(rz, rd) / crossprod(rz)); den <- drop(crossprod(dh, rd)); b <- drop(crossprod(dh, ry)) / den
      e <- drop(ry - b * rd); g <- dh * e; L <- h - 1
      v0 <- sum(g^2); for (a in seq_len(L)) v0 <- v0 + 2 * (L + 1 - a) / (L + 1) * sum(g[(a + 1):nh] * g[1:(nh - a)])
      res[h, ] <- c(b, sqrt(v0) / abs(den)) }
  }
  res
}
