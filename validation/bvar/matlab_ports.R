# R ports of the deterministic parts of naffe15/BVAR_ (bvartools/*.m). One function per .m file.
# Conventions follow MATLAB: y is T x ny, Phi is (ny*lags+nx) x ny with lags first then constant.
suppressMessages(library(MASS))
BVAR_DIR <- "/Users/VASILK04/Documents/research/external/BVAR_"
EXDIR <- file.path(BVAR_DIR, "examples", "BVAR tutorial")
source("/Users/VASILK04/Documents/research/varm/validation/bvar/readmat5.R")
# varm is loaded with pkgload::load_all() from VARM_SRC (a scratch copy of the repo, so that compiling never touches the
# shared working tree; the repo's own load_all raced with other agents' compiles). Default: the repo itself.
loadvarm <- function() {
  src <- Sys.getenv("VARM_SRC", "/Users/VASILK04/Documents/research/varm")
  suppressMessages(pkgload::load_all(src, quiet = TRUE)); invisible(TRUE)
}
rep_ <- function(name, ok, ...) {
  cat(sprintf("[%s] %s  %s\n", if (is.na(ok)) "NOT-TESTED" else if (ok) "PASS" else "FAIL", name, paste(..., collapse = " ")))
  invisible(ok)
}
mx <- function(a, b) max(abs(a - b))

# ---- rfvar3.m (Sims) --------------------------------------------------------------------------
rfvar3 <- function(ydata, lags, xdata, breaks, lambda, mu, plr = NULL, ybar = NULL) {
  nvar <- ncol(ydata); T <- nrow(ydata); nx <- ncol(xdata)
  breaks <- c(0, breaks, T); smpl <- integer(0)
  for (nb in seq_len(length(breaks) - 1)) if (breaks[nb] + lags + 1 <= breaks[nb + 1]) smpl <- c(smpl, seq(breaks[nb] + lags + 1, breaks[nb + 1]))
  X <- matrix(0, length(smpl), nvar * lags)
  for (is in seq_along(smpl)) for (l in 1:lags) X[is, (l - 1) * nvar + 1:nvar] <- ydata[smpl[is] - l, ]
  X <- cbind(X, xdata[smpl, , drop = FALSE]); y <- ydata[smpl, , drop = FALSE]
  if (lambda != 0 || mu > 0) {
    if (is.null(ybar)) ybar <- colMeans(ydata[1:lags, , drop = FALSE]); xbar <- if (nx) colMeans(xdata[1:lags, , drop = FALSE]) else numeric(0)
    if (lambda != 0) {
      xdum <- if (lambda > 0) lambda * c(rep(ybar, lags), xbar) else { lambda <- -lambda; lambda * c(rep(ybar, lags), 0 * xbar) }
      y <- rbind(y, lambda * ybar); X <- rbind(X, xdum)
    }
    if (mu > 0) {
      xd <- cbind(do.call(cbind, rep(list(diag(ybar, nvar)), lags)), matrix(0, nvar, nx)) * mu
      X <- rbind(X, xd); y <- rbind(y, mu * diag(ybar, nvar))
    }
  }
  s <- svd(X); di <- 1 / s$d
  B <- (s$v * rep(di, each = nvar * lags + nx)) %*% t(s$u) %*% y
  xxi <- s$v * rep(di, each = nvar * lags + nx); xxi <- xxi %*% t(xxi)
  list(B = B, u = y - X %*% B, xxi = xxi, y = y, X = X)
}
# ---- varprior.m ---------------------------------------------------------------------------------
varprior <- function(nv, nx, lags, tight, decay, sig, w) {
  K3 <- lags * nv
  yd <- array(0, c(lags + 1, nv, K3)); xd <- array(0, c(lags + 1, nx, K3))
  for (j in 1:nv) for (il in 1:lags) {
    k <- il + lags * (j - 1)
    yd[il + 1, j, k] <- il^decay * sig[j]
    if (il == 1) yd[1, j, k] <- sig[j]
  }
  yd <- tight * yd; yd <- yd[(lags + 1):1, , , drop = FALSE]
  breaks <- (lags + 1) * seq_len(K3); lbreak <- breaks[K3]
  if (w > 0) for (i in 1:w) {
    y2 <- array(0, c(lags + 1, nv, nv)); for (j in 1:nv) y2[lags + 1, j, j] <- sig[j]
    yd <- array(c(yd, y2), c(lags + 1, nv, dim(yd)[3] + nv)); xd <- array(0, c(lags + 1, nx, dim(yd)[3]))
    breaks <- c(breaks, (lags + 1) * seq_len(nv) + lbreak); lbreak <- breaks[length(breaks)]
  }
  d3 <- dim(yd)[3]
  Y <- matrix(0, (lags + 1) * d3, nv); for (k in 1:d3) Y[(k - 1) * (lags + 1) + 1:(lags + 1), ] <- yd[, , k]
  Xd <- matrix(0, (lags + 1) * d3, nx)
  list(ydum = Y, xdum = Xd, breaks = breaks[-length(breaks)])
}
# ---- matrictint.m -------------------------------------------------------------------------------
matrictint <- function(S, df, XXi) {
  k <- nrow(XXi); ny <- nrow(S); cx <- chol(XXi); cs <- chol(S)
  w1 <- 0.5 * k * ny * log(2 * pi) + ny * sum(log(diag(cx)))
  w2 <- -df * sum(log(diag(cs))) + 0.5 * df * ny * log(2) + ny * (ny - 1) * 0.25 * log(pi) +
    sum(lgamma(0.5 * (df + (0:(1 - ny)))))
  w1 + w2
}
# ---- bvar_ml.m (Minnesota dummy prior, presample = 0, first_obs = lags+1, constant) --------------
# hyper = c(tau, decay, lambda, mu, omega). Returns list(logml, post, prior).
bvar_ml <- function(hyper, y, lags, constant = TRUE, ybar_shift = 0) {   # ybar_shift = p mimics varm's ybar rows
  ny <- ncol(y); nx <- as.integer(constant)
  tau <- hyper[1]; decay <- hyper[2]; lambda <- hyper[3]; mu <- hyper[4]; omega <- hyper[5]
  sig <- apply(y[1:lags, , drop = FALSE], 2, sd)            # std(y(1:lags,:)) -- MATLAB std uses n-1
  vp <- varprior(ny, nx, lags, tau, decay, sig, round(omega))
  ydata <- y; T <- nrow(ydata); xdata <- matrix(1, T, nx)
  yb <- if (ybar_shift) colMeans(y[ybar_shift + 1:lags, , drop = FALSE]) else NULL
  post <- rfvar3(rbind(ydata, vp$ydum), lags, rbind(xdata, vp$xdum), c(T, T + vp$breaks), lambda, mu, ybar = yb)
  Tp <- lags
  prior <- rfvar3(rbind(ydata[1:Tp, , drop = FALSE], vp$ydum), lags, rbind(xdata[1:Tp, , drop = FALSE], vp$xdum),
                  c(Tp, Tp + vp$breaks), lambda, mu, ybar = yb)
  pdf <- nrow(prior$u) - ny * lags - nx; qdf <- nrow(post$u) - ny * lags - nx
  pi_ <- matrictint(crossprod(post$u), qdf, post$xxi); pr_ <- matrictint(crossprod(prior$u), pdf, prior$xxi)
  list(logml = pi_ - pr_ - 0.5 * ny * (qdf - pdf) * log(2 * pi), post = c(post, list(df = qdf)), prior = c(prior, list(df = pdf)), sig = sig)
}
# ---- iresponse.m -------------------------------------------------------------------------------
# alpha: (ny*lags) x ny. Returns array [var, horizon(1..hor), shock]. unit=1: one std; unit=0: Q*inv(diag(diag(Q))) (1 pct)
iresponse <- function(alpha, Sigma, hor, Omega, unit = 1) {
  N <- nrow(Sigma); lags <- nrow(alpha) %/% N
  Q <- t(chol(Sigma)); if (unit == 0) Q <- Q %*% solve(diag(diag(Q), N))
  F <- rbind(t(alpha[1:(N * lags), ]), if (lags > 1) cbind(diag(N * (lags - 1)), matrix(0, N * (lags - 1), N)))
  G <- diag(N * lags)[, 1:N, drop = FALSE]
  ir <- array(0, c(N, hor, ncol(Omega))); Fk <- diag(N * lags)
  for (k in 1:hor) { ir[, k, ] <- t(G) %*% (Fk %*% G %*% Q %*% Omega); Fk <- F %*% Fk }
  ir
}
# ---- fevd.m -------------------------------------------------------------------------------------
fevd_m <- function(hor, Phi, Sigma, Omega = diag(nrow(Sigma))) {
  N <- nrow(Sigma); lags <- nrow(Phi) %/% N
  F <- rbind(t(Phi[1:(N * lags), ]), if (lags > 1) cbind(diag(N * (lags - 1)), matrix(0, N * (lags - 1), N)))
  G <- diag(N * lags)[, 1:N, drop = FALSE]; A <- t(chol(Sigma)); Kp <- G %*% A %*% Omega
  Fp <- diag(N * lags); tot <- 0; V <- matrix(0, N, N)
  for (hh in 1:hor) {
    for (s in 1:N) { tt <- Fp %*% Kp[, s, drop = FALSE]; V[, s] <- V[, s] + diag(tt %*% t(tt))[1:N] }
    Fp <- F %*% Fp
  }
  V / rowSums(V) * 100
}
# ---- iresponse_longrun.m ------------------------------------------------------------------------
iresponse_longrun <- function(alpha, Sigma, hor, lags) {
  N <- nrow(Sigma)
  F <- rbind(t(alpha[1:(N * lags), ]), if (lags > 1) cbind(diag(N * (lags - 1)), matrix(0, N * (lags - 1), N)))
  C1 <- solve(diag(N * lags) - F)[1:N, 1:N]
  PSI1 <- t(chol(C1 %*% Sigma %*% t(C1)))
  Q <- solve(C1) %*% PSI1
  Q[1, 1] <- Q[1, 1] * sign(Q[1, 1])               # NB: only element (1,1), as in the .m file
  G <- diag(N * lags)[, 1:N, drop = FALSE]; ir <- array(0, c(N, hor, N)); Fk <- diag(N * lags)
  for (k in 1:hor) { ir[, k, ] <- t(G) %*% (Fk %*% G %*% Q); Fk <- F %*% Fk }
  list(ir = ir, Q = Q)
}
# ---- iresponse_proxy.m (point estimate, one instrument, constant in both stages) ----------------
iresponse_proxy <- function(res, proxy, Phi, Sigma, p, irhor, T_m_end = 0) {
  T <- nrow(res); n <- ncol(res); Tm <- NROW(proxy); k <- 1
  ix <- (T - Tm - T_m_end + 1):(T - T_m_end); Z <- cbind(1, proxy)
  Phib <- solve(crossprod(Z), crossprod(Z, res[ix, ])); uhat1 <- Z %*% Phib[, 1]
  b21ib11 <- t(solve(crossprod(cbind(1, uhat1)), crossprod(cbind(1, uhat1), res[ix, -1, drop = FALSE]))[2, , drop = FALSE])
  Sig11 <- Sigma[1, 1]; Sig21 <- Sigma[-1, 1, drop = FALSE]; Sig22 <- Sigma[-1, -1, drop = FALSE]
  ZZp <- b21ib11 %*% Sig11 %*% t(b21ib11) - (Sig21 %*% t(b21ib11) + b21ib11 %*% t(Sig21)) + Sig22
  b12b12p <- t(Sig21 - b21ib11 %*% Sig11) %*% solve(ZZp, Sig21 - b21ib11 %*% Sig11)
  b11 <- sqrt(Sig11 - b12b12p)
  b1 <- c(b11, b21ib11 * drop(b11))
  irs <- matrix(0, p + irhor, n); irs[p + 1, ] <- b1
  for (jj in 2:irhor) {
    lv <- t(irs[(p + jj - 1):jj, , drop = FALSE]); irs[p + jj, ] <- as.vector(lv) %*% Phi[1:(p * n), ]
  }
  list(irs = irs[p + 1:irhor, ], b1 = b1, Phib = Phib, uhat1 = uhat1)
}
# ---- forecasts.m (no-shock path) ----------------------------------------------------------------
forecast_noshock <- function(init, Phi, fhor, xdata) {   # init: lags x ny (oldest first)
  ny <- ncol(init); lags <- nrow(init); out <- matrix(0, fhor, ny); ld <- init
  for (t in 1:fhor) {
    X <- c(as.vector(t(ld[nrow(ld):1, , drop = FALSE])), xdata[t, ]); y <- X %*% Phi
    ld[1:(lags - 1), ] <- ld[-1, ]; ld[lags, ] <- y; out[t, ] <- y
  }
  out
}
# ---- YXB_.m / ols (rfvar3 with no dummies, lambda=mu=0) ----------------------------------------
yxb <- function(y, lags, constant = TRUE) {
  nobs <- nrow(y) - lags; X <- matrix(0, nobs, ncol(y) * lags)
  for (i in 1:lags) X[, (i - 1) * ncol(y) + 1:ncol(y)] <- y[(lags - i + 1):(lags + nobs - i), ]
  if (constant) X <- cbind(X, 1)
  list(Y = y[lags + 1:nobs, , drop = FALSE], X = X)
}
# ---- histdecomp.m (Omega = I or given; yo = first lags obs, constant only) ----------------------
histdecomp <- function(u, alpha, Sigma, data, lags, Omega = diag(ncol(Sigma))) {
  N <- ncol(Sigma); A <- t(chol(Sigma)); ierror <- u %*% solve(t(Omega) %*% t(A))
  F <- rbind(t(alpha[1:(N * lags), ]), if (lags > 1) cbind(diag(N * (lags - 1)), matrix(0, N * (lags - 1), N)))
  G <- diag(N * lags)[, 1:N, drop = FALSE]
  C <- c(alpha[N * lags + 1, ], rep(0, N * (lags - 1)))
  yo <- as.vector(t(data[lags:1, ]))                # == reshape(flip(...)) = XX(1,1:N*lags)
  Tu <- nrow(u); Fp <- vector("list", Tu + 1); Fp[[1]] <- diag(N * lags); for (i in 1:Tu) Fp[[i + 1]] <- F %*% Fp[[i]]
  Kp <- G %*% A %*% Omega; E <- array(0, c(Tu, N, N)); Bd <- matrix(0, Tu, N)
  for (t in 1:Tu) {
    Aa <- 0; for (tau in 1:t) Aa <- Aa + Fp[[tau]] %*% C
    Bd[t, ] <- (Aa + Fp[[t + 1]] %*% yo)[1:N]
    for (s in 1:N) { D <- 0; for (tau in 1:t) D <- D + Fp[[t - tau + 1]] %*% Kp[, s, drop = FALSE] * ierror[tau, s]; E[t, , s] <- D[1:N] }
  }
  list(E = E, B = Bd, yhat = Bd + apply(E, c(1, 2), sum))
}
