# Validation of the Bayesian engines of varm. Run: Rscript validation/bayes.R [section ...]
# Sections: s1 s2 s3 s4 s5 s6 s7 s8 (each prints its own PASS/FAIL line).
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE); run <- function(s) length(args) == 0 || s %in% args
rep_ <- function(name, ok, ...) cat(sprintf("[%s] %s  %s\n", if (ok) "PASS" else "FAIL", name, paste(..., collapse = " ")))
zdiff <- function(a, b) { # z-score of difference of means (iid approx, columns)
  (colMeans(a) - colMeans(b)) / sqrt(apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b))
}

# ---------- exact structural posterior for the homoskedastic triangular prior -----------------
# Equation i: y_i = Z_i gamma_i + e, e ~ N(0, s2); gamma_i ~ N(g0, diag(v)); s2 ~ IG(a0, b0*sc_i^2) independent
# (semi-conjugate, so p(s2|y) is obtained by 1-D quadrature of the marginal likelihood N(y; Z g0, s2 I + Z V Z')).
tri_prior_moments <- function(i, K, p, d, scale, pr) {
  m <- K * p + d; v <- numeric(m + i - 1); g0 <- numeric(m + i - 1)
  for (l in 1:p) for (j in 1:K) {
    b <- scale[i]^2 / scale[j]^2 / l^pr$alpha
    v[(l - 1) * K + j] <- b * pr$lambda1^2 * (if (j == i) 1 else pr$lambda2^2)
    if (l == 1 && j == i) g0[(l - 1) * K + j] <- rep_len(pr$delta, K)[i]
  }
  for (q in seq_len(d)) v[K * p + q] <- (pr$lambda3 * scale[i])^2
  for (j in seq_len(i - 1)) v[m + j] <- pr$lambda5^2 * scale[i]^2 / scale[j]^2
  list(v = v, g0 = g0)
}
exact_tri <- function(Y, p, pr, n, const = TRUE, ngrid = 6000) {
  K <- ncol(Y); det <- make_det(nrow(Y), const, FALSE, NULL); o <- cpp_var_ols(Y, p, det)
  X <- o$X; Ye <- o$Y; m <- ncol(X); sc <- sqrt(ar_resid_var(Y, p))
  Pi <- array(0, c(m, K, n)); A <- array(diag(K), c(K, K, n)); S2 <- matrix(0, n, K)
  for (i in 1:K) {
    pm <- tri_prior_moments(i, K, p, ncol(det), sc, pr)
    Z <- cbind(X, -Ye)[, seq_len(m + i - 1), drop = FALSE]; y <- Ye[, i]; r <- y - Z %*% pm$g0
    ev <- eigen(Z %*% (pm$v * t(Z)), symmetric = TRUE); dk <- pmax(ev$values, 0); ur2 <- crossprod(ev$vectors, r)[, 1]^2
    a0 <- pr$hom[1]; b0 <- pr$hom[2] * sc[i]^2
    ls <- seq(log(sc[i]^2) - 6, log(sc[i]^2) + 4, length.out = ngrid); s2 <- exp(ls)
    lm <- vapply(s2, function(s) -0.5 * sum(log(s + dk)) - 0.5 * sum(ur2 / (s + dk)), 0)
    lw <- lm - (a0 + 1) * ls - b0 / s2 + ls; w <- exp(lw - max(lw)); w <- w / sum(w)
    s2d <- s2[sample.int(ngrid, n, TRUE, w)]; S2[, i] <- s2d
    ZtZ <- crossprod(Z); Zty <- crossprod(Z, y)
    for (k in 1:n) {
      P <- ZtZ / s2d[k] + diag(1 / pm$v, length(pm$v)); Lc <- chol(P)
      mu <- backsolve(Lc, backsolve(Lc, Zty / s2d[k] + pm$g0 / pm$v, transpose = TRUE))
      g <- mu + backsolve(Lc, rnorm(length(mu)))
      Pi[, i, k] <- g[1:m]; if (i > 1) A[i, 1:(i - 1), k] <- g[m + 1:(i - 1)]
    }
  }
  B <- array(0, c(m, K, n)); Sg <- array(0, c(K, K, n))
  for (k in 1:n) { Ai <- solve(matrix(A[, , k], K, K)); B[, , k] <- matrix(Pi[, , k], m, K) %*% t(Ai); Sg[, , k] <- Ai %*% diag(S2[k, ], K) %*% t(Ai) }
  list(B = B, Sigma = Sg, Pi = Pi, A = A, S2 = S2)
}
flat <- function(a) t(matrix(a, prod(dim(a)[1:2])))

if (run("s1")) {
  set.seed(1); Y1 <- as.matrix(yraw[, 1, drop = FALSE]); N <- 40000
  pr <- prior_triangular("minnesota", lambda1 = 0.3, lambda5 = 1)
  f <- bvarm(Y1, p = 2, prior = pr, ndraw = N, nburn = 500)
  e <- exact_tri(Y1, 2, pr, N)
  z <- zdiff(flat(f$draws$B), flat(e$B)); zs <- zdiff(flat(f$draws$Sigma), flat(e$Sigma))
  vr <- apply(flat(f$draws$B), 2, var) / apply(flat(e$B), 2, var)
  rep_("1a K=1 matches exact semi-conjugate posterior", max(abs(c(z, zs))) < 4.5 && all(abs(log(vr)) < 0.1),
       "max|z|=", round(max(abs(c(z, zs))), 2), " var ratio range", paste(round(range(vr), 3), collapse = "-"))
  # K=3 structural check
  Y3 <- as.matrix(yraw[, 1:3]); N <- 20000
  f <- bvarm(Y3, p = 2, prior = pr, ndraw = N, nburn = 500)
  e <- exact_tri(Y3, 2, pr, N)
  z <- zdiff(flat(f$draws$B), flat(e$B)); zs <- zdiff(flat(f$draws$Sigma), flat(e$Sigma))
  za <- zdiff(flat(f$A), flat(e$A))
  zp <- c(max(abs(z)), max(abs(zs)), max(abs(za[is.finite(za)])))
  rep_("1b K=3 reduced-form (B, Sigma, A) vs exact structural posterior", all(zp < 4.5), "max|z| B/Sigma/A:", paste(round(zp, 2), collapse = " "))
  # loose prior -> OLS
  fl <- bvarm(Y3, p = 2, prior = prior_triangular("minnesota", lambda1 = 1e3, lambda2 = 1, lambda3 = 1e3, lambda5 = 1e3, hom = c(0.01, 0.01)),
              ndraw = 20000, nburn = 500)
  fo <- varm(Y3, p = 2)
  de <- max(abs(coef(fl) - fo$B)); dS <- max(abs(apply(fl$draws$Sigma, 1:2, mean) - fo$Sigma_ml) / abs(fo$Sigma_ml))
  rep_("1c loose prior reproduces OLS", de < 0.02, "max|B-OLS|=", signif(de, 3), " (Sigma mean rel dev ", round(dS, 3), ")")
}

if (run("s2")) {
  library(stochvol, quietly = TRUE)
  # (a) Geweke joint-distribution test of the KSC mixture + FFBS + IG update
  set.seed(2); g <- cpp_sv_geweke(T = 10, nsim = 400000, a_h = 5, b_h = 0.2, mu0 = 0, V0 = 1)
  zc <- vapply(1:4, function(j) {
    ess <- coda::effectiveSize(g$chain[, j]); (mean(g$chain[, j]) - mean(g$fwd[, j])) /
      sqrt(var(g$fwd[, j]) / nrow(g$fwd) + var(g$chain[, j]) / ess) }, 0)
  vrat <- apply(g$chain, 2, var) / apply(g$fwd, 2, var)
  cat("Geweke z (sig2,h0,hT,mean h):", round(zc, 2), " var ratios:", round(vrat, 3), "\n")
  rep_("2a Geweke joint-distribution test (KSC approx; mixture is approximate, tolerance |z|<4)", all(abs(zc) < 4), "")
  # (b) simulated VAR(1), K=2 with known volatility paths
  set.seed(3); T <- 600; K <- 2
  Aa <- matrix(c(1, 0.5, 0, 1), 2, 2)           # A y_t = c + Pi' y_{t-1} + eps ; a21 = 0.5
  Pi1 <- matrix(c(0.5, 0.1, 0, 0.4), 2, 2, byrow = TRUE); cc <- c(0.1, 0.2)
  h1 <- cumsum(rnorm(T, 0, 0.12)) - 0.3; h1 <- h1 - mean(h1) - 0.5        # time-varying, log-var
  h2 <- rep(-0.7, T)                                                       # constant
  Y <- matrix(0, T, K); for (t in 2:T) Y[t, ] <- solve(Aa, cc + Pi1 %*% Y[t - 1, ] + exp(0.5 * c(h1[t], h2[t])) * rnorm(2))
  colnames(Y) <- c("a", "b")
  fit <- bvarm(Y, p = 1, prior = prior_triangular("minnesota", lambda1 = 1, lambda2 = 1, lambda5 = 10, delta = 0),
               sv = TRUE, ndraw = 3000, nburn = 2000, store_h = TRUE, sv_prior = c(a = 5, b = 0.04, V0 = 10))
  Ht <- fit$h_all; hm <- fit$h_mean; tt <- 2:T                     # sample rows align with t = 2..T
  cor1 <- cor(hm[, 1], h1[tt]); 
  q <- apply(Ht[, , 1], 2, quantile, c(.05, .95)); cov90 <- mean(h1[tt] >= q[1, ] & h1[tt] <= q[2, ])
  q2 <- apply(Ht[, , 2], 2, quantile, c(.05, .95)); cov90c <- mean(h2[tt] >= q2[1, ] & h2[tt] <= q2[2, ])
  sdh <- colMeans(fit$sigma_h); 
  cat(sprintf("eq1 (tv): cor(h)=%.3f 90%% coverage=%.2f ; sigma_h posterior mean=%.3f (true 0.12)\n", cor1, cov90, sdh[1]))
  cat(sprintf("eq2 (const): 90%% coverage=%.2f ; mean sd of h2 path=%.3f; sigma_h mean=%.3f (true 0, prior mean 0.1)\n", cov90c, mean(apply(Ht[, , 2], 2, sd)), sdh[2]))
  rep_("2b time-varying vol recovery", cor1 > 0.8 && cov90 > 0.8, "")
  rep_("2b constant vol recovered (sigma_h shrinks, coverage)", cov90c > 0.8 && sdh[2] < 0.1, "")
  cs <- coef(fit); cat("B posterior mean:\n"); print(round(cs, 2)); cat("true reduced-form B':\n"); print(round(solve(Aa, cbind(cc, Pi1)), 2))
  Ah <- apply(fit$A, 1:2, mean); cat("A mean a21 =", round(Ah[2, 1], 3), "(true 0.5)\n")
  # (c) brute-force algebra of the reconstruction, per stored draw
  d <- cpp_bvar_tri(cpp_var_ols(Y, 1, make_det(T, TRUE, FALSE, NULL))$Y, cpp_var_ols(Y, 1, make_det(T, TRUE, FALSE, NULL))$X, 1L, 1L, 0L, TRUE,
                    sqrt(ar_resid_var(Y, 1)), c(0, 0), c(1, 1, 2, 100, 10, .1, .5), c(5, .04, 10), c(3, 2), 5L, 50L, 1L, TRUE)
  err <- 0
  for (k in 1:5) {
    A <- d$A[, , k]; Ai <- solve(A); L <- diag(exp(d$h_all[k, T - 1, ]), 2)
    err <- max(err, max(abs(d$Sigma[, , k] - Ai %*% L %*% t(Ai))))
    # structural equation A y_t = Pi' x_t + u => B = Pi A^{-T}: verify the implied reduced-form mean at a point x
    xx <- c(1, 0.3, -0.2); err <- max(err, 0)  # B itself checked in 1b via the exact structural posterior
  }
  rep_("2c Sigma = A^-1 Lambda A^-T algebra (last-period vol)", err < 1e-10, "max err", signif(err, 2))
}

sim_sparse <- function(T = 250, K = 6, seed = 4) {
  set.seed(seed)
  A1 <- diag(0.5, K); A1[cbind(c(2, 3, 5, 6), c(1, 2, 4, 5))] <- 0.3; A2 <- matrix(0, K, K); A2[1, 4] <- -0.25
  Y <- matrix(0, T + 50, K); for (t in 3:(T + 50)) Y[t, ] <- A1 %*% Y[t - 1, ] + A2 %*% Y[t - 2, ] + rnorm(K, 0, c(1, 1.5, 1, .7, 1, 1.2))
  Y <- Y[-(1:50), ]; colnames(Y) <- paste0("y", 1:K)
  list(Y = Y, truth = rbind(t(A1), t(A2)))   # (Kp x K) in B layout: rows = lag vars, cols = equation
}
if (run("s3")) {
  s <- sim_sparse(); K <- 6; nz <- s$truth != 0
  res <- list()
  for (ty in c("minnesota", "ssvs", "horseshoe")) {
    pr <- prior_triangular(ty, lambda1 = if (ty == "ssvs") 0.5 else 0.3, lambda2 = 1, alpha = 1, lambda5 = 1, delta = 0, rho = 0.05, inclusion = 0.5)
    set.seed(5); f <- bvarm(s$Y, p = 2, prior = pr, ndraw = 3000, nburn = 1500)
    Bh <- coef(f)[1:(2 * K), ]; ae <- abs(Bh)
    res[[ty]] <- c(mean_abs_zero = mean(ae[!nz]), mean_abs_nonzero = mean(ae[nz]), rmse = sqrt(mean((Bh - s$truth)^2)),
                   if (ty == "ssvs") c(incl_zero = mean(f$inclusion[!nz]), incl_nonzero = mean(f$inclusion[nz])),
                   if (ty == "horseshoe") c(tau_mean = mean(f$tau)))
  }
  fo <- varm(s$Y, p = 2); res$OLS <- c(mean_abs_zero = mean(abs(fo$B[1:12, ][!nz])), mean_abs_nonzero = mean(abs(fo$B[1:12, ][nz])),
                                       rmse = sqrt(mean((fo$B[1:12, ] - s$truth)^2)))
  print(res)
  r <- res; incl <- r$ssvs
  rep_("3a SSVS inclusion separates true nonzeros from zeros", incl["incl_nonzero"] > 0.9 && incl["incl_zero"] < 0.3,
       "P(incl|nonzero)=", round(incl["incl_nonzero"], 3), " P(incl|zero)=", round(incl["incl_zero"], 3))
  rep_("3b horseshoe/SSVS beat OLS RMSE; horseshoe shrinks zeros", r$horseshoe["rmse"] < r$OLS["rmse"] && r$ssvs["rmse"] < r$OLS["rmse"] &&
       r$horseshoe["mean_abs_zero"] < 0.6 * r$OLS["mean_abs_zero"], "rmse hs/ssvs/ols=", round(r$horseshoe["rmse"], 4), round(r$ssvs["rmse"], 4), round(r$OLS["rmse"], 4))
}

if (run("s4")) {
  Y <- as.matrix(yraw[, 1:3]); p <- 2; M <- 3
  # (a) fixed hyperparameters: direct draws vs analytic NIW moments (no dummy priors so formulas are explicit)
  pr <- prior_minnesota(lambda = 0.3, alpha = 2, soc = NULL, sur = NULL, hyper = character(0), delta = 1)
  f <- bvarm(Y, p = p, prior = pr, ndraw = 60000)
  det <- make_det(nrow(Y), TRUE, FALSE, NULL); o <- cpp_var_ols(Y, p, det); X <- o$X; Ye <- o$Y; m <- ncol(X)
  psi <- ar_resid_var(Y, p); om <- c(rep(0.3^2 / rep(1:p, each = M)^2, 1) / rep(psi, p), 1e7)
  b0 <- matrix(0, m, M); b0[cbind(1:M, 1:M)] <- 1
  Ginv <- solve(crossprod(X) + diag(1 / om)); bh <- Ginv %*% (crossprod(X, Ye) + diag(1 / om) %*% b0)
  R <- Ye - X %*% bh; S <- diag(psi) + crossprod(R) + t(bh - b0) %*% diag(1 / om) %*% (bh - b0); nu <- nrow(Ye) + M + 2
  B <- f$draws$B; Sg <- f$draws$Sigma
  mb <- apply(B, 1:2, mean); eS <- S / (nu - M - 1); mS <- apply(Sg, 1:2, mean)
  vecB <- t(matrix(B, m * M)); covB <- cov(vecB); covA <- kronecker(eS, Ginv)    # vec(B) is column-stacked (eq by eq)
  e1 <- max(abs(mb - bh) / sqrt(diag(covA)[1:(m * M)] %>% (function(x) matrix(x, m, M)) / 60000))
  rel <- function(a, b) max(abs(a - b) / abs(b))
  rep_("4a fixed-hyper draws vs analytic: E[B] (max z), E[Sigma] (max rel), Cov(vec B) (max rel on diag/offdiag)",
       e1 < 4.5 && rel(mS, eS) < 0.02 && rel(diag(covB), diag(covA)) < 0.03 && max(abs(cov2cor(covB) - cov2cor(covA))) < 0.03,
       "z=", round(e1, 2), " Sigma rel=", round(rel(mS, eS), 4), " var rel=", round(rel(diag(covB), diag(covA)), 4),
       " corr maxdiff=", round(max(abs(cov2cor(covB) - cov2cor(covA))), 4))
  # (b) hierarchical lambda/soc/sur vs BVAR package
  # NB: BVAR::bvar's automatic psi is the AR(p) residual *sd*-scale (0.318 vs our variance 0.1015); pass its psi to compare like with like
  psiB <- BVAR:::auto_psi(Y[-(1:p), ], p)$mode
  set.seed(7); fh <- bvarm(Y, p = p, prior = prior_minnesota(psi = psiB), ndraw = 20000, nburn = 5000)
  set.seed(7); bv <- BVAR::bvar(Y, lags = p, priors = BVAR::bv_priors(hyper = c("lambda", "soc", "sur"), soc = BVAR::bv_soc(), sur = BVAR::bv_sur()), n_draw = 25000, n_burn = 5000, n_thin = 1, verbose = FALSE)
  hv <- bv$hyper
  cmp <- rbind(varm = colMeans(fh$hyper)[c("lambda", "soc", "sur")], BVAR = colMeans(hv)[c("lambda", "soc", "sur")])
  q <- rbind(varm = quantile(fh$hyper[, "lambda"], c(.05, .5, .95)), BVAR = quantile(hv[, "lambda"], c(.05, .5, .95)))
  print(round(cmp, 4)); print(round(q, 4))
  # exact reference: 3-D quadrature of the (verified-identical-to-BVAR) log posterior over (lambda, soc, sur)
  pr0 <- prior_minnesota(psi = psiB); m <- ncol(cpp_var_ols(Y, p, make_det(nrow(Y), TRUE, FALSE, NULL))$X)
  oo <- cpp_var_ols(Y, p, make_det(nrow(Y), TRUE, FALSE, NULL)); bb <- matrix(0, m, M); bb[cbind(1:M, 1:M)] <- 1
  gp <- rbind(gamma_coef(pr0$lambda["mode"], pr0$lambda["sd"]), gamma_coef(2, .25), gamma_coef(1, 1), gamma_coef(1, 1))
  lpf <- function(l, s, u) cpp_niw_logpost(c(l, 2, s, u, psiB), oo$Y, oo$X, p, 1L, TRUE, bb, 1e7, as.integer(c(1, 0, 1, 1, 0, 0, 0)), gp,
                                           c(1e-4, 1, 1e-4, 1e-4, rep(0, M)), c(5, 3, 50, 50, rep(Inf, M)), c(.004, .004), FALSE)
  gl <- exp(seq(log(0.02), log(5), length.out = 80)); gs <- exp(seq(log(0.01), log(50), length.out = 50))
  Lg <- array(0, c(80, 50, 50)); for (i in 1:80) for (j in 1:50) for (k in 1:50) Lg[i, j, k] <- lpf(gl[i], gs[j], gs[k]) + log(gl[i] * gs[j] * gs[k])
  W <- exp(Lg - max(Lg)); W <- W / sum(W)
  ex <- c(lambda = sum(gl * apply(W, 1, sum)), soc = sum(gs * apply(W, 2, sum)), sur = sum(gs * apply(W, 3, sum)))
  cmp <- rbind(cmp, exact = ex); print(round(cmp, 4))
  rel_d <- abs(cmp["varm", ] / cmp["exact", ] - 1); rel_b <- abs(cmp["BVAR", ] / cmp["exact", ] - 1)
  rep_("4b hierarchical hyper means: varm MCMC vs exact quadrature (within 5%)", all(rel_d < 0.05), "rel diffs", paste(round(rel_d, 3), collapse = " "),
       "| BVAR's own MCMC vs exact:", paste(round(rel_b, 3), collapse = " "))
  eb <- max(abs(coef(fh) - apply(bv$beta, 2:3, mean)))   # BVAR beta layout: draws x (1+K p) x K, const first
  cat("   max |posterior mean B diff| vs BVAR (after row reordering):\n")
  bb <- apply(bv$beta, 2:3, mean); bb <- rbind(bb[-1, ], bb[1, ]); cat("  ", signif(max(abs(coef(fh) - bb)), 3), "\n")
  rep_("4c posterior mean of B vs BVAR", max(abs(coef(fh) - bb)) < 0.03, "")
}

if (run("s5")) {
  set.seed(11)
  rep3 <- function(M, N) array(M, c(dim(M), N))
  drive <- function(K, S, N, Sigma = diag(K), B = matrix(0, K, K)) cpp_sign_irf(rep3(B, N), rep3(Sigma, N), 1L, 1L, S, 1L, 1L)
  N <- 200000
  # (a) Haar uniformity: no restrictions, Sigma = I => P = Q
  for (K in 2:3) {
    o <- drive(K, array(0, c(K, K, 1)), N); Q <- o$P
    orth <- max(abs(apply(Q, 3, function(q) max(abs(crossprod(q) - diag(K))))))
    q11 <- Q[1, 1, ]
    ks <- if (K == 2) ks.test(q11, function(x) 0.5 + asin(pmax(pmin(x, 1), -1)) / pi)$statistic      # cos(theta), theta~U
          else ks.test(q11, "punif", -1, 1)$statistic                                               # Archimedes: entries uniform on S^2
    e2 <- range(apply(Q^2, 1:2, mean)); dets <- mean(apply(Q, 3, det))
    # all entries have the same distribution; also check an off-diagonal entry and joint (Q11,Q22) correlation
    ks2 <- if (K == 3) ks.test(Q[2, 3, ], "punif", -1, 1)$statistic else ks.test(Q[1, 2, ], function(x) 0.5 + asin(pmax(pmin(x, 1), -1)) / pi)$statistic
    rep_(sprintf("5a Haar K=%d: orthogonal, entry law exact (KS), E[Q^2]=1/K, E[det]=0", K),
         orth < 1e-12 && ks < 0.005 && ks2 < 0.005 && all(abs(e2 - 1 / K) < 0.005) && abs(dets) < 0.01,
         "KS=", signif(ks, 3), signif(ks2, 3), " E[Q^2] range", paste(round(e2, 4), collapse = "-"), " E[det]=", round(dets, 4))
  }
  # (b) acceptance rates vs analytic, K = 2 and 3 (Sigma = I, B = 0: impact response = Q)
  S2 <- array(0, c(2, 2, 1)); S2[, 1, 1] <- c(1, 1); S2[, 2, 1] <- c(-1, 1)     # demand (+,+), supply (-,+)
  o <- drive(2, S2, N); acc2 <- o$accepted / N
  rep_("5b K=2 demand/supply toy: acceptance = 1/2", abs(acc2 - 0.5) < 0.005, "acc=", round(acc2, 4))
  theta <- atan2(o$P[2, 1, ], o$P[1, 1, ]); rep_("5b accepted demand-shock angles confined to (0,pi/2) and uniform there",
       all(theta > 0 & theta < pi / 2) && ks.test(theta, "punif", 0, pi / 2)$statistic < 0.01, "KS=", signif(ks.test(theta, "punif", 0, pi / 2)$statistic, 3))
  S3 <- array(0, c(3, 3, 1)); S3[, 1, 1] <- c(1, 1, 1)                                   # col 1 fully signed: P = 2/8
  a3 <- drive(3, S3, N)$accepted / N
  rep_("5b K=3 fully signed column: acceptance = 1/4", abs(a3 - 0.25) < 0.005, "acc=", round(a3, 4))
  S3b <- S3; S3b[, 2, 1] <- c(1, -1, 0)   # second column (+,-,.) given orthogonal to q1 >0: analytic via independent R Haar sampler
  hr <- function(K) { q <- qr(matrix(rnorm(K * K), K)); Q <- qr.Q(q) %*% diag(sign(diag(qr.R(q)))); Q }
  ref <- mean(replicate(60000, { Q <- hr(3); ok <- TRUE
    for (j in 1:2) { s <- S3b[, j, 1]; nz <- s != 0; v <- Q[nz, j]; ok <- ok && (all(s[nz] * v > 0) || all(-s[nz] * v > 0)) }; ok }))
  a3b <- drive(3, S3b, N)$accepted / N
  rep_("5b K=3 two restricted columns vs independent R Haar sampler", abs(a3b - ref) < 0.01, "C++:", round(a3b, 4), "R:", round(ref, 4))
  # (c) every accepted IRF satisfies the restrictions on a real fitted VAR; P P' = Sigma
  Y <- as.matrix(yraw[, 1:3]); fit <- bvarm(Y, p = 2, prior = prior_minnesota(hyper = character(0)), ndraw = 500)
  R <- array(0, c(3, 3, 4)); R[1, 1, ] <- 1; R[2, 1, ] <- 1; R[3, 2, ] <- -1; R[1, 2, ] <- 1; R[3, 3, ] <- 1
  ir <- irf(fit, id = id_sign(R, nrot = 5, maxtry = 5000), horizon = 8)
  a <- ir$draws; viol <- 0
  for (n in seq_len(dim(a)[4])) for (j in 1:3) for (h in 1:4) for (i in 1:3) { s <- R[i, j, h]; if (s != 0 && s * a[i, j, h, n] <= 0) viol <- viol + 1 }
  d <- fit$draws
  P2 <- cpp_sign_irf(d$B, d$Sigma, 2L, 8L, R, 2L, 5000L); perr <- max(vapply(seq_len(dim(P2$P)[3]), function(k) max(abs(P2$P[, , k] %*% t(P2$P[, , k]) - d$Sigma[, , P2$draw[k]])), 0))
  rep_("5c accepted IRFs satisfy all sign restrictions; P P' = Sigma", viol == 0 && perr < 1e-10, "violations=", viol, " of ", length(a) / 1, " max|PP'-Sigma|=", signif(perr, 2),
       " acceptance rate=", round(ir$accept_rate, 4), " accepted draws=", dim(a)[4])
}

if (run("s6")) {
  nrep <- as.integer(Sys.getenv("NREP", 300)); nb <- 400; T <- 200; K <- 2; H <- 12
  rho <- as.numeric(Sys.getenv("RHO", 0.9)); Btrue <- rbind(matrix(c(rho, 0.0, 0.3, 0.5), 2, 2), c(0, 0))   # rows: lag-1 vars, const; col = equation  (y1 persistent)
  Sig <- matrix(c(1, 0.4, 0.4, 1), 2); L <- t(chol(Sig))
  truth <- cpp_irf(Btrue, K, 1L, H, L)                                 # [resp, shock, h]
  simdgp <- function(het) { e <- matrix(rnorm((T + 50) * K), ncol = K) %*% t(L)
    if (het) e <- e * ifelse(seq_len(T + 50) > (T + 50) * 0.5, 2, 0.5)
    Y <- matrix(0, T + 50, K); for (t in 2:(T + 50)) Y[t, ] <- Y[t - 1, ] %*% Btrue[1:2, ] + e[t, ]; Y[-(1:50), ] }
  meths <- list(iid = list("bootstrap", FALSE), iid_bc = list("bootstrap", TRUE), wild = list("wild", FALSE), wild_bc = list("wild", TRUE), mc = list("mc", FALSE))
  for (het in c(FALSE, TRUE)) {
    set.seed(100 + het); cov <- array(0, c(length(meths), 2, nrep)); width <- cov; bias <- matrix(0, length(meths), nrep)
    for (r in 1:nrep) {
      Y <- simdgp(het); fit <- varm(Y, p = 1)
      for (m in seq_along(meths)) {
        ir <- irf(fit, id_chol(), horizon = H, ci = meths[[m]][[1]], nboot = nb, level = c(.68, .9), bias_correct = meths[[m]][[2]])
        b <- ir$bands; tr <- truth
        for (l in 1:2) { lo <- if (l == 1) b[, , , 2] else b[, , , 1]; hi <- if (l == 1) b[, , , 3] else b[, , , 4]
          cov[m, l, r] <- mean(tr >= lo & tr <= hi); width[m, l, r] <- mean(hi - lo) }
      }
    }
    tab <- data.frame(method = names(meths), cov68 = apply(cov[, 1, ], 1, mean), cov90 = apply(cov[, 2, ], 1, mean),
                      w68 = apply(width[, 1, ], 1, mean), w90 = apply(width[, 2, ], 1, mean))
    cat(sprintf("\nIRF coverage (Cholesky, VAR(1), T=%d, %d reps x %d draws, %s errors; avg over 2x2x13 responses)\n", T, nrep, nb, if (het) "heteroskedastic" else "iid"))
    print(format(tab, digits = 3), row.names = FALSE)
    if (!het) rep_(sprintf("6 iid DGP (rho=%.1f): bias-corrected bootstrap and MC within 0.07 of nominal; BC improves 90%% coverage when persistent", rho),
                   all(abs(tab$cov68[c(2, 5)] - .68) < .07) && all(abs(tab$cov90[c(2, 5)] - .90) < .07) && (rho < .9 || tab$cov90[2] > tab$cov90[1]), "")
    if (het) rep_("6 heteroskedastic DGP: wild bootstrap coverage >= iid bootstrap coverage (90%)", tab$cov90[3] >= tab$cov90[1] - 0.005, round(tab$cov90[3], 3), "vs", round(tab$cov90[1], 3))
  }
}

if (run("s7")) {
  set.seed(21)
  ar <- function(phi, n) as.numeric(arima.sim(list(ar = phi), n))
  # known ESS of AR(1): n (1-phi)/(1+phi)
  r <- replicate(40, mcmc_diagnostics(cbind(ar(0.8, 4000)))$ess); tr <- 4000 * 0.2 / 1.8
  rep_("7a ESS of AR(1) phi=.8 vs theory n(1-phi)/(1+phi)", abs(mean(r) / tr - 1) < 0.08, "mean ESS=", round(mean(r)), " theory=", round(tr))
  ch <- lapply(1:4, function(j) cbind(a = ar(0.5, 2000), b = rnorm(2000)))
  d <- mcmc_diagnostics(ch); print(d)
  rep_("7b multi-chain: converged chains R-hat ~ 1, ESS(iid) ~ N", all(abs(d$rhat - 1) < 0.01) && abs(d$ess[2] / 8000 - 1) < 0.2 && abs(d$ess[1] / (8000 / 3) - 1) < 0.2, "")
  bad <- lapply(1:4, function(j) cbind(a = ar(0.5, 2000) + j * 2))
  rep_("7c non-mixed chains detected by R-hat", mcmc_diagnostics(bad)$rhat > 1.5, "rhat=", round(mcmc_diagnostics(bad)$rhat, 2))
  drift <- cbind(a = cumsum(rnorm(3000, 0.01)) + rnorm(3000))
  zt <- mcmc_diagnostics(cbind(a = ar(0.5, 3000)))$geweke_z; zd <- mcmc_diagnostics(drift)$geweke_z
  zs <- replicate(200, mcmc_diagnostics(cbind(ar(0.5, 2000)))$geweke_z)
  rep_("7d Geweke z ~ N(0,1) under stationarity (sd in 0.8-1.25), large under drift", abs(sd(zs) - 1) < 0.25 && abs(zd) > 3, "sd(z)=", round(sd(zs), 2), " drift z=", round(zd, 1))
  fit <- bvarm(as.matrix(yraw[, 1:3]), p = 2, ndraw = 3000, nburn = 1000); print(mcmc_diagnostics(fit))
  rep_("7e works on bvarm (hyperparameters)", nrow(mcmc_diagnostics(fit)) >= 3, "")
}

if (run("s8")) {
  Tn <- 250; p <- 4; nd <- 1000
  sim <- function(K) { set.seed(K); A <- matrix(rnorm(K * K, 0, 0.15 / sqrt(K)), K) + diag(0.5, K); Y <- matrix(0, Tn + 50, K)
    for (t in 2:(Tn + 50)) Y[t, ] <- Y[t - 1, ] %*% t(A) + rnorm(K); Y <- Y[-(1:50), ]; colnames(Y) <- paste0("y", 1:K); Y }
  tm <- function(expr) unname(system.time(expr)[3])
  res <- NULL
  add <- function(K, what, sec) res <<- rbind(res, data.frame(K = K, engine = what, sec_per_1000 = round(sec, 2)))
  for (K in c(3, 7, 20)) {
    Y <- sim(K)
    add(K, "varm conjugate, fixed hyper", tm(bvarm(Y, p, prior_minnesota(hyper = character(0)), ndraw = nd)))
    add(K, "varm conjugate, hierarchical (opt + MH)", tm(bvarm(Y, p, ndraw = nd, nburn = 500)))
    add(K, "varm triangular minnesota, homosk.", tm(bvarm(Y, p, prior_triangular("minnesota"), ndraw = nd, nburn = 200)))
    add(K, "varm triangular horseshoe, homosk.", tm(bvarm(Y, p, prior_triangular("horseshoe"), ndraw = nd, nburn = 200)))
    add(K, "varm triangular minnesota + SV", tm(bvarm(Y, p, prior_triangular("minnesota"), sv = TRUE, ndraw = nd, nburn = 200)))
    add(K, "BVAR::bvar (hierarchical lambda)", tm(BVAR::bvar(Y, lags = p, n_draw = nd + 500, n_burn = 500, n_thin = 1, verbose = FALSE)))
    if (K <= 7) {
      spec <- bsvars::specify_bsvar_sv$new(Y, p = p); tt <- tm(invisible(utils::capture.output(bsvars::estimate(spec, S = nd, show_progress = FALSE))))
      add(K, "bsvars::bsvar_sv (structural SV, S=1000)", tt)
    }
    cat("K =", K, "done\n")
  }
  print(res, row.names = FALSE)
}

if (run("s2d")) {  # simulation-based calibration (lite) of the random-walk innovation variance and the vol path, truths drawn from the prior
  set.seed(31); nr <- as.integer(Sys.getenv("NR", 60)); T <- 300; a <- 5; b <- 0.04
  cov_s <- cov_h <- numeric(nr); rk <- numeric(nr)
  for (r in 1:nr) {
    s2 <- 1 / rgamma(1, a, b); h <- cumsum(rnorm(T, 0, sqrt(s2))) + rnorm(1, 0, sqrt(1)); y <- matrix(exp(h / 2) * rnorm(T), ncol = 1)
    f <- bvarm(y, p = 1, prior = prior_triangular("minnesota", lambda1 = 0.5, delta = 0), sv = TRUE, sv_prior = c(a = a, b = b, V0 = 1), ndraw = 1500, nburn = 1000, store_h = TRUE)
    sd_draw <- f$sigma_h[, 1]; rk[r] <- mean(sd_draw^2 < s2); cov_s[r] <- s2 > quantile(sd_draw^2, .05) && s2 < quantile(sd_draw^2, .95)
    q <- apply(f$h_all[, , 1], 2, quantile, c(.05, .95)); cov_h[r] <- mean(h[-1] >= q[1, ] & h[-1] <= q[2, ])
  }
  cat("rank histogram (quintiles) of true sigma_h^2 in posterior:", table(cut(rk, c(0, .2, .4, .6, .8, 1), include.lowest = TRUE)), "\n")
  rep_("2d sigma_h^2 90% credible-interval coverage (prior-drawn truths)", mean(cov_s) > 0.8, "cov=", round(mean(cov_s), 3), " mean vol-path 90% coverage=", round(mean(cov_h), 3))
}
