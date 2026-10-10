# Default (non-informative) reduced-form posterior: BVAR_ bvar_.m with no Minnesota prior vs varm draws(fit, "mc")
source("/Users/VASILK04/Documents/research/varm/validation/bvar/matlab_ports.R"); loadvarm()
set.seed(20261011)
GK <- read_mat5(file.path(EXDIR, "DataGK.mat")); y <- cbind(logip = GK$logip, logcpi = GK$logcpi, gs1 = GK$gs1, ebp = GK$ebp); lags <- 2; K <- 4
o <- yxb(y, lags); Te <- nrow(o$Y); k <- ncol(o$X); Phi <- qr.coef(qr(o$X, tol = 1e-12), o$Y); u <- o$Y - o$X %*% Phi; S <- crossprod(u); XXi <- chol2inv(qr.R(qr(o$X)))
# BVAR_: flat = 1 => posterior.df = Tu - ny*lags - nx - (ny+1);  Sigma ~ rand_inverse_wishart(ny, df, chol(inv(S)))  (IW(S, df): mean S/(df - ny - 1));  Phi | Sigma ~ N(PhiHat, Sigma (x) XXi)
dfB <- Te - k - (K + 1)
dfV <- Te - k       # varm cpp_var_mc: nu = Te - m  (IW(U'U, nu))
cat("   Te =", Te, " k =", k, " K =", K, ": BVAR_ df =", dfB, " varm nu =", dfV, "\n")
fit <- varm(y, lags); N <- 40000
dv <- draws(fit, "mc", N)
draw_B <- function(df, N) { Sc <- chol(solve(S)); Ss <- array(0, c(K, K, N)); for (n in 1:N) { X <- matrix(rnorm(df * K), df, K) %*% Sc; s <- svd(X); Ss[, , n] <- (s$v * rep(1 / s$d^2, each = K)) %*% t(s$v) }; Ss }
SB <- draw_B(dfB, N)
mS <- function(a) apply(a, 1:2, mean)
exp_B <- S / (dfB - K - 1); exp_V <- S / (dfV - K - 1)
zB <- (mS(SB) - exp_B)[cbind(1:4, 1:4)] / sqrt(apply(SB, c(1, 2), var)[cbind(1:4, 1:4)] / N); zV <- (mS(dv$Sigma) - exp_V)[cbind(1:4, 1:4)] / sqrt(apply(dv$Sigma, c(1, 2), var)[cbind(1:4, 1:4)] / N)
rep_("D1 BVAR_ sampler: E[Sigma] = S/(df-K-1) with df = Te-k-(K+1)", max(abs(zB)) < 4.5, "max|z| diag=", round(max(abs(zB)), 2))
rep_("D2 varm mc sampler: E[Sigma] = U'U/(nu-K-1) with nu = Te-m (Jeffreys posterior)", max(abs(zV)) < 4.5, "max|z| diag=", round(max(abs(zV)), 2))
cat("   ratio E_varm[Sigma]/E_BVAR_[Sigma] (diag) =", signif(diag(exp_V / exp_B), 5), " (=", signif((dfB - K - 1) / (dfV - K - 1), 5), ")\n")
rep_("D3 default-prior Sigma posteriors coincide (BVAR_ df vs varm nu)", dfB == dfV, "BVAR_ df=", dfB, " varm nu=", dfV, " -> BVAR_ posterior is IW(S, T-k-K-1) (flat prior p(Sigma)=1, although labelled Jeffreys); varm is IW(S, T-k) (Jeffreys)")
# B | Sigma draws have the same structure: Cov(vec B | Sigma) = Sigma (x) XXi in both
zb <- sapply(1:3, function(j) { i <- sample(1:ncol(dv$B[, , 1]), 1); NA })
Bm <- t(matrix(dv$B, k * K)); Vth <- kronecker(mS(dv$Sigma), XXi)    # unconditional cov = E[Sigma] (x) XXi
rep_("D4 varm mc: Cov(vec B) = E[Sigma] (x) (X'X)^-1 (diagonal, 3% tolerance)", max(abs(diag(cov(Bm)) / diag(Vth) - 1)) < 0.03, "max rel dev=", signif(max(abs(diag(cov(Bm)) / diag(Vth) - 1)), 3))
# ---- IC: BVAR_ IC.m (AIC/HQIC/BIC with S = dof-corrected Sigma_ols, K = regressors per equation, own sample per lag) vs varm var_select (common sample, Lutkepohl form) ----
IC <- function(S, E, T, Kx) { N <- nrow(S); llf <- -(T * N / 2) * (1 + log(2 * pi)) - T / 2 * log(det(S)) - 0.5 * sum(diag(MASS::ginv(S) %*% crossprod(E)))
  c(AIC = -2 * llf / T + 2 * Kx / T, HQIC = -2 * llf / T + 2 * Kx * log(log(T)) / T, BIC = -2 * llf / T + Kx * log(T) / T) }
icB <- sapply(1:6, function(p) { oo <- yxb(y, p); Ph <- qr.coef(qr(oo$X, tol = 1e-12), oo$Y); e <- oo$Y - oo$X %*% Ph; nk <- ncol(oo$X); IC(crossprod(e) / (nrow(e) - nk), e, nrow(e), nk) })
vs <- var_select(y, 6)
cat("   BVAR_ IC.m argmin lags (AIC, HQIC, BIC):", apply(icB, 1, which.min), " ; varm var_select (AIC, HQ, SC):", vs$selection[c("AIC", "HQ", "SC")], "\n")
rep_("I1 information-criterion lag choice agrees (BVAR_ IC.m vs varm var_select, 1..6 lags)", all(apply(icB, 1, which.min) == vs$selection[c("AIC", "HQ", "SC")]),
     "BVAR_:", paste(apply(icB, 1, which.min), collapse = ","), " varm:", paste(vs$selection[c("AIC", "HQ", "SC")], collapse = ","), " (BVAR_ penalty uses K regressors/equation, not K*N + N(N+1)/2; own sample per lag)")
vv <- vars::VARselect(y, lag.max = 6, type = "const")$selection
rep_("I2 varm var_select agrees with vars::VARselect (AIC, HQ, SC)", all(vv[c("AIC(n)", "HQ(n)", "SC(n)")] == vs$selection[c("AIC", "HQ", "SC")]), "vars:", paste(vv[1:3], collapse = ","), " varm:", paste(vs$selection[c("AIC", "HQ", "SC")], collapse = ","))
