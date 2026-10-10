# Heavier validation of R/lp-extra.R. Run from the package root: Rscript validation/lp-extra.R
pkgload::load_all(".", quiet = TRUE)
df <- as.data.frame(yraw); H <- 8
cat("== 1. lp_state(timing = 'lagged') vs lpirfs::lp_nl ==\n")
# Conventions: state F(z_{t-1}) with z standardised (HP cycle or raw), common constant, lag-1 coefficients x Cholesky
# column of the linear VAR. lpirfs bands are sum_i d_i (b_i +- se_i) (NOT delta-method), so bands are
# compared through the per-coefficient HAC se returned in $raw.
for (cfg in list(list(hp = 1, sh = 0, tr = 0, gam = 3), list(hp = 0, sh = 1, tr = 1, gam = 2))) {
  zz <- if (cfg$hp == 0) as.numeric(scale(yraw$infl)) else yraw$infl
  ref <- lpirfs::lp_nl(df, lags_endog_lin = 2, lags_endog_nl = 2, trend = cfg$tr, shock_type = cfg$sh, confint = 1, hor = H,
                       switching = zz, use_hp = cfg$hp, lambda = 1600, gamma = cfg$gam, use_nw = TRUE, num_cores = 1)
  for (sj in 1:3) {
    me <- lp_state(yraw, sj, "infl", H, lags = 2, timing = "lagged", lags_lin = 2, shock_type = cfg$sh, gamma = cfg$gam,
                   hp_lambda = if (cfg$hp) 1600, state_const = FALSE, trend = cfg$tr == 1)
    d <- me$contrasts[1:3, "high"]; mxse <- 0
    for (k in 1:3) for (h in 1:H) {
      V <- me$raw$vcov[(h - 1) * 6 + 1:6, (h - 1) * 6 + 1:6, k]; s <- sqrt(diag(V))
      mxse <- max(mxse, abs(sum(d * s[1:3]) - (ref$irf_s1_up[k, h + 1, sj] - ref$irf_s1_mean[k, h + 1, sj])),
                  abs(sum(d * s[4:6]) - (ref$irf_s2_up[k, h + 1, sj] - ref$irf_s2_mean[k, h + 1, sj])))
    }
    cat(sprintf("hp=%d shock_type=%d trend=%d shock %d: max|IRF diff| high %.1e low %.1e; max|band half-width diff| %.1e\n", cfg$hp, cfg$sh, cfg$tr, sj,
                max(abs(me$estimate[, , "high"] - t(ref$irf_s1_mean[, , sj]))), max(abs(me$estimate[, , "low"] - t(ref$irf_s2_mean[, , sj]))), mxse))
  }
}

cat("\n== 2. sup-t coverage and smooth LP MSE (simulated VAR(1), 300 reps, T=200) ==\n")
A <- matrix(c(0.6, 0.2, 0, 0.5), 2, 2); P <- t(chol(matrix(c(1, .3, .3, 1), 2)))
Ht <- 12; true <- sapply(0:Ht, function(h) { M <- diag(2); for (i in seq_len(h)) M <- A %*% M; (M %*% P[, 1] / P[1, 1])[2] })
set.seed(42); R <- 300; cov_pw <- cov_st <- 0; mse_raw <- mse_sm <- 0; crit <- numeric(R); cnt <- 0
for (r in seq_len(R)) {
  e <- matrix(rnorm(2 * 250), 250) %*% t(P); Y <- matrix(0, 250, 2)
  for (t in 2:250) Y[t, ] <- A %*% Y[t - 1, ] + e[t, ]
  Y <- Y[51:250, ]; colnames(Y) <- c("x", "y")
  f <- lp_joint(Y, "x", Ht, lags = 2, responses = "y", level = 0.9)
  b <- lp_supt(f, 0.9, nsim = 2000, seed = r); crit[r] <- attr(b, "crit")[[1]]
  cov_pw <- cov_pw + all(true >= b$pw_lower & true <= b$pw_upper); cov_st <- cov_st + all(true >= b$lower & true <= b$upper)
  sm <- lp_smooth(Y, "x", "y", Ht, lags = 2)
  mse_raw <- mse_raw + mean((f$estimate[, 1] - true)^2); mse_sm <- mse_sm + mean((sm$estimate[, 1] - true)^2)
}
cat(sprintf("90%% simultaneous coverage of the whole path: pointwise bands %.2f, sup-t bands %.2f (nominal 0.90); mean sup-t crit %.2f\n", cov_pw / R, cov_st / R, mean(crit)))
cat(sprintf("MSE of the response path: raw LP %.4f, smooth LP %.4f (ratio %.2f)\n", mse_raw / R, mse_sm / R, mse_sm / mse_raw))

cat("\n== 2b. same design, lag-augmented (HC) errors: pointwise vs sup-t coverage ==\n")
set.seed(42); cv <- matrix(0, R, Ht + 1); st <- numeric(R)
for (r in seq_len(R)) {
  e <- matrix(rnorm(2 * 250), 250) %*% t(P); Y <- matrix(0, 250, 2)
  for (t in 2:250) Y[t, ] <- A %*% Y[t - 1, ] + e[t, ]
  Y <- Y[51:250, ]; colnames(Y) <- c("x", "y")
  b <- lp_supt(lp_joint(Y, "x", Ht, lags = 1, responses = "y", lag_augment = TRUE), 0.9, nsim = 2000, seed = r)
  cv[r, ] <- true >= b$pw_lower & true <= b$pw_upper; st[r] <- all(true >= b$lower & true <= b$upper)
}
cat(sprintf("mean pointwise coverage %.2f; whole-path sup-t coverage %.2f (HAC/HC finite-sample undercoverage carries over)\n", mean(cv), mean(st)))

cat("\n== 3. LP vs VAR (Plagborg-Moller & Wolf): same lags, recursive ID, T=5000 ==\n")
set.seed(7); e <- matrix(rnorm(2 * 5050), 5050) %*% t(P); Y <- matrix(0, 5050, 2)
for (t in 2:5050) Y[t, ] <- A %*% Y[t - 1, ] + e[t, ]; Y <- Y[-(1:50), ]; colnames(Y) <- c("x", "y")
cmp <- compare_lp_var(Y, "x", 12, lags = 1)
cat(sprintf("max |LP - VAR| over horizons 0..12 = %.3f; max |VAR - truth| = %.3f; max |LP - truth| = %.3f\n",
            max(abs(cmp$diff)), max(abs(cmp$var[cmp$response == "y"] - true)), max(abs(cmp$lp[cmp$response == "y"] - true))))
cmp <- compare_lp_var(yraw, "unemp", 12, lags = 2); print(cmp)

cat("\n== 4. lp_panel speed (N=200, T=60, two-way FE, H=8) ==\n")
set.seed(3); N <- 200; Tn <- 60; d <- expand.grid(id = 1:N, t = 1:Tn); d$s <- rnorm(nrow(d))
d$y <- 0.5 * d$s + rnorm(N)[d$id] + rnorm(Tn)[d$t] + rnorm(nrow(d)); d <- d[order(d$id, d$t), ]
print(system.time(r <- lp_panel(d, "y", "s", "id", "t", horizon = 8, lags = 2)))
print(round(cbind(est = r$estimate[, 1], se = r$se[, 1]), 3)[1:3, ])
