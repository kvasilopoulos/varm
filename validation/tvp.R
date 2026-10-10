# Validation of tvpvar(): run with the package installed (or pkgload::load_all) -- takes a few minutes.
#  1. constant-coefficient limit (tv_beta = FALSE, and kQ -> 0) vs OLS on the estimation sample
#  2. DK-smoother posterior sd (kQ -> 0) vs exact GLS posterior sd of the constant model
#  3. coverage of true paths on simulated TVP-SV DGPs
#  4. timing T = 215, K = 3, 10000 draws
suppressMessages(if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) pkgload::load_all(".") else library(varm))
set.seed(42)

# ---- 1 & 2 ----------------------------------------------------------------
p <- 2; tau <- 40
fit_ols <- varm(yraw[(p + tau + 1 - p):nrow(yraw), ], p = p)           # estimation sample incl. p lags
pr0 <- tvp_prior(tau = tau, kQ = 1e-6, kS = 1e-6, kW = 1e-6, vscale = 1e4)   # diffuse initial states: isolates the likelihood
f_c  <- tvpvar(yraw, p, ndraw = 2000, nburn = 1000, sv = FALSE, tv_beta = FALSE, tv_A = FALSE, prior = pr0)
f_q0 <- tvpvar(yraw, p, ndraw = 2000, nburn = 1000, sv = FALSE, tv_A = FALSE, prior = pr0)
B_ols <- fit_ols$B
cat("max |posterior mean (constant) - OLS|     :", format(max(abs(coef(f_c, date = "last") - B_ols)), digits = 3), "\n")
cat("max |posterior mean (kQ = 1e-6) - OLS|    :", format(max(abs(coef(f_q0, date = "last") - B_ols)), digits = 3), "\n")
cat("max |time variation of mean path (kQ->0)| :", format(max(apply(f_q0$beta_mean, 1:2, function(v) diff(range(v)))), digits = 3), "\n")
sdc <- apply(f_c$beta[, , 1, ], 1:2, sd); sdq <- apply(f_q0$beta[, , 87, ], 1:2, sd)
cat("posterior sd, DK (kQ->0) / exact GLS: range ", format(range(sdq / sdc), digits = 3), "\n")
f_sv0 <- tvpvar(yraw, p, ndraw = 2000, nburn = 1000, sv = TRUE, prior = pr0)
cat("with SV on, kQ,kS,kW -> 0: max |mean - OLS| :", format(max(abs(coef(f_sv0, date = "last") - B_ols)), digits = 3),
    " (volatility sd range ", format(range(f_sv0$w_sd), digits = 2), ")\n")
# BVAR comparison (weak Minnesota / flat) on the same sample
bv <- bvarm(yraw[(tau + 1):nrow(yraw), ], p = p, prior = prior_minnesota(lambda = c(5, 1, 1e-4, 50), soc = NULL, sur = NULL, hyper = character(0)),
            ndraw = 2000)
cat("max |BVAR(loose) mean - OLS| :", format(max(abs(apply(bv$draws$B, 1:2, mean) - B_ols)), digits = 3), "\n")

# ---- 3a exact-prior calibration (simulate from the prior used by the sampler; Q, S, W pinned by huge prior df)
cal <- function(rep = 20, T = 150, tv_beta = TRUE, tv_A = TRUE, sv = TRUE) {
  K <- 2; m <- 3; nb <- 6
  b0 <- c(0.5, 0.1, 0, 0.1, 0.5, 0); P0 <- diag(0.01, nb); Qv <- 1e-4; Sv <- 1e-3; Wv <- 0.01; a0 <- 0.3
  res <- matrix(NA, rep, 3)
  for (r in 1:rep) {
    b <- b0 + 0.1 * rnorm(nb); a <- a0 + 0.2 * rnorm(1); h <- log(0.25) + sqrt(0.5) * rnorm(2)
    Y <- matrix(0, T, K); X <- matrix(0, T, m); y <- c(0, 0); TB <- matrix(0, nb, T); TA <- numeric(T); TH <- matrix(0, 2, T)
    for (t in 1:T) {
      if (tv_beta) b <- b + sqrt(Qv) * rnorm(nb); if (tv_A) a <- a + sqrt(Sv) * rnorm(1); if (sv) h <- h + sqrt(Wv) * rnorm(2)
      x <- c(y, 1); X[t, ] <- x
      A <- matrix(c(1, 0, a, 1), 2, 2, byrow = TRUE)
      y <- as.numeric(x %*% matrix(b, m, K) + solve(A, exp(h / 2) * rnorm(2))); Y[t, ] <- y
      TB[, t] <- b; TA[t] <- a; TH[, t] <- h
    }
    o <- varm:::cpp_tvp(Y, X, 1L, b0, P0, a0, list(NULL, matrix(0.04)), diag(Qv * 1e6, nb), 1e6, diag(Qv, nb),
         list(NULL, matrix(Sv * 1e6)), c(0, 1e6), list(NULL, matrix(Sv)), rep(log(0.25), 2), 0.5, c(1e5, 1e5 * Wv),
         c(3, 2), c(sv, tv_beta, tv_A, 0, 0, 1), 400L, 200L, 1L, 400L)
    cv <- function(a, tr) { q <- apply(a, 1:2, quantile, probs = c(.05, .95)); mean(tr >= q[1,,] & tr <= q[2,,]) }
    res[r, ] <- c(cv(array(o$beta, c(nb, T, 400)), TB), cv(array(o$alpha, c(1, T, 400)), matrix(TA, 1)), cv(o$h, TH))
  }
  colMeans(res)
}
set.seed(3); cc <- cal(); cat(sprintf("calibration, nominal 0.90: beta %.3f, alpha %.3f, log-vol %.3f (20 reps, T = 150)\n", cc[1], cc[2], cc[3]))

# ---- 3b coverage with the data-driven training-sample prior (prior != DGP: shows prior-tightness effect) -------------------------------------------------------------
sim_tvp <- function(T = 250, K = 3, p = 1, qb = 0.012, qa = 0.03, qh = 0.08) {
  m <- K * p + 1
  B0 <- rbind(diag(0.6, K), c(0.2, 0.1, 0.1)); nb <- m * K
  beta <- matrix(0, nb, T); alpha <- matrix(0, K * (K - 1) / 2, T); h <- matrix(0, K, T)
  b <- as.numeric(B0); a <- c(0.3, -0.2, 0.2); hh <- log(c(0.3, 0.2, 0.25)^2)
  Y <- matrix(0, T, K); y <- rep(0, K)
  for (t in 1:T) {
    b <- b + qb * rnorm(nb); a <- a + qa * rnorm(length(a)); hh <- hh + qh * rnorm(K)
    A <- diag(K); A[2, 1] <- a[1]; A[3, 1] <- a[2]; A[3, 2] <- a[3]
    Bt <- matrix(b, m, K)
    y <- as.numeric(c(y, 1) %*% Bt + solve(A, exp(hh / 2) * rnorm(K)))
    Y[t, ] <- y; beta[, t] <- b; alpha[, t] <- a; h[, t] <- hh
  }
  list(Y = Y, beta = beta, alpha = alpha, h = h)
}
for (kQ in c(0.01, 0.05)) {
nrep <- 8; cov_b <- cov_h <- cov_a <- rmse_b <- rmse_const <- numeric(0)
for (r in 1:nrep) {
  d <- sim_tvp(); if (max(abs(d$Y)) > 50) next
  colnames(d$Y) <- c("a", "b", "c")
  f <- tvpvar(d$Y, p = 1, ndraw = 1500, nburn = 1000, prior = tvp_prior(tau = 40, kQ = kQ))
  est <- 41:250 + 0   # p = 1: effective rows 2..250, est rows after tau -> original rows 42..250
  est <- (1 + 40 + 1):250
  tb <- d$beta[, est]; bq <- apply(f$beta, 1:3, quantile, probs = c(0.05, 0.95))
  tbm <- array(tb, dim(f$beta)[1:3])
  cov_b <- c(cov_b, mean(tbm >= bq[1, , , ] & tbm <= bq[2, , , ]))
  hq <- apply(f$h, 1:2, quantile, probs = c(0.05, 0.95)); th <- d$h[, est]
  cov_h <- c(cov_h, mean(th >= hq[1, , ] & th <= hq[2, , ]))
  aq <- apply(f$alpha, 1:2, quantile, probs = c(0.05, 0.95)); ta <- d$alpha[, est]
  cov_a <- c(cov_a, mean(ta >= aq[1, , ] & ta <= aq[2, , ]))
  rmse_b <- c(rmse_b, sqrt(mean((f$beta_mean - tbm)^2)))
  rmse_const <- c(rmse_const, sqrt(mean((rowMeans(tb) - as.numeric(f$beta_mean[, , 1]) * 0 - tb)^2)))
}
cat(sprintf("kQ = %g: coverage of 90%% bands over %d reps: beta %.3f, log-variance h %.3f, A_t %.3f\n", kQ, length(cov_b), mean(cov_b), mean(cov_h), mean(cov_a)))
cat(sprintf("RMSE of posterior-mean beta path %.4f vs RMSE of best constant (true mean) %.4f\n", mean(rmse_b), mean(rmse_const)))
}

# ---- 4 timing ---------------------------------------------------------------
tm <- system.time(ff <- tvpvar(yraw, p = 2, ndraw = 10000, nburn = 2000))
cat("T = 215, K = 3, p = 2, SV + TV-A + full Q, 10000 draws (+2000 burn), 1000 stored paths: ", tm[["elapsed"]], " s\n")
tm <- system.time(tvpvar(yraw, p = 4, ndraw = 10000, nburn = 2000))
cat("same with p = 4 (state dim 39): ", tm[["elapsed"]], " s\n")
