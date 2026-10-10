# Identification through heteroskedasticity: varm id_hetero() vs BVAR_ iresponse_heterosked.m.  Simulated VAR(1) with known structural impact matrix P.
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm()
set.seed(20261011)
K <- 3; T <- 4000; P0 <- matrix(c(1, .4, -.3, 0, 1, .5, .3, -.2, 1), K) %*% diag(c(1, 1.2, .8))   # true impact matrix (columns = shocks)
lam <- c(4, 1.6, .5)                                                                                   # regime-2 / regime-1 variance ratios (distinct)
A1 <- matrix(c(.5, .1, 0, 0, .4, .1, .1, 0, .3), K)
reg <- rep(0:1, each = T / 2); eps <- matrix(rnorm(T * K), T) * outer(rep(1, T), ifelse(rep(1, K) > 0, 1, 1)); eps[reg == 1, ] <- sweep(eps[reg == 1, ], 2, sqrt(lam), "*")
y <- matrix(0, T, K); for (t in 2:T) y[t, ] <- y[t - 1, ] %*% t(A1) + t(P0 %*% eps[t, ])
colnames(y) <- c("a", "b", "c")
fit <- varm(y, 1); id <- id_hetero(fit, reg)
Pv <- varm:::structural_irf(array(fit$B, c(dim(fit$B), 1)), array(fit$Sigma, c(K, K, 1)), K, 1L, 0L, id, colnames(y))$irf[, , 1, 1]
match_cols <- function(P, P0) { # align by max |cosine| and sign; return max relative error of aligned columns
  cs <- abs(cor(P, P0)); perm <- apply(cs, 1, which.max); if (anyDuplicated(perm)) return(Inf)
  Pa <- P0[, perm]; Pa <- Pa %*% diag(sign(colSums(P * Pa))); max(abs(P - Pa)) / max(abs(P0)) }
ev <- match_cols(Pv, P0)
rep_("HT1 varm id_hetero recovers the true impact matrix (T=4000, up to order/sign)", ev < 0.1, "max rel. error=", signif(ev, 3))
rep_("HT2 varm id_hetero: lambda (variance ratios) recovered", max(abs(sort(attr(id, "lambda")) - sort(lam))) < 0.3, "est=", paste(round(attr(id, "lambda"), 3), collapse = ","), " true=", paste(lam, collapse = ","))
# exact population check: with exact covariances S1 = P P', S2 = P L P' the varm formula returns P exactly
S1 <- P0 %*% t(P0); S2 <- P0 %*% diag(lam) %*% t(P0); L1 <- t(chol(S1)); e <- eigen(solve(L1, t(solve(L1, S2))), symmetric = TRUE); Pp <- L1 %*% e$vectors
rep_("HT3 varm formula on exact population covariances returns P (error ~ 1e-14)", match_cols(Pp, P0) < 1e-10, "err=", signif(match_cols(Pp, P0), 3))
# ---- BVAR_ iresponse_heterosked.m on the same data ----
u <- fit$U; d <- rep(reg, 1)[-1]; S0 <- cov(u[d == 0, ]); S1h <- cov(u[d == 1, ]); X <- solve(S0, S1h); sv <- svd(X)
Wb <- diag(sqrt(sv$d)) %*% sv$v                          # iresponse(Phi, S, hor, V): impact = chol(S)' * Omega with S = diag(singular values), Omega = V
cat("   BVAR_ impact (chol(diag(svals)) * V):\n"); print(round(Wb, 3)); cat("   true P0:\n"); print(round(P0, 3))
rep_("HT4 BVAR_ impact matrix reproduces the sample covariance (W W' = Sigma_regime1)", mx(Wb %*% t(Wb), S0) / max(abs(S0)) < 0.05, "rel err=", signif(mx(Wb %*% t(Wb), S0) / max(abs(S0)), 3), " (W W' = diag(svals) by construction)")
rep_("HT5 BVAR_ impact matrix recovers the true P (up to order/sign)", match_cols(Wb, P0) < 0.1, "max rel. error=", signif(match_cols(Wb, P0), 3))
# what a correct use of the same svd would give: impact = chol(S0) * V' is NOT right either; the identified matrix is chol(S0) * eigvec(chol(S0)^-1 S1 chol(S0)^-T)
rep_("HT6 population check of BVAR_ formula: exact S0=PP', S1=P L P' -> BVAR_ impact", match_cols({ s <- svd(solve(S1, S2)); diag(sqrt(s$d)) %*% s$v }, P0) < 1e-6, "(expected FAIL if the .m code is wrong)")
