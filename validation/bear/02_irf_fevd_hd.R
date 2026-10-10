# BEAR vs varm: Cholesky IRF, FEVD, HD on the OLS VAR (BEAR irfsim.m / fevd.m / hdecomp.m replicated in R)
source("validation/bear/common.R")
p <- 4; d <- load_bear_data(p); Y <- d$Y; X <- d$X; K <- 3; H <- 13  # BEAR identificationHorizon=12 => 13 points incl. impact
f <- varm(Y, p, const = TRUE, exog = X); o <- bear_olsvar(Y, X, TRUE, p)
rn <- rownames(f$B); idx <- c(1:(K * p), which(rn == "const"), which(rn == "Oil"))
D <- t(chol(o$Sigma))                           # BEAR: Sigma = D D', D lower (Cholesky identifier)
ir <- bear_irfsim(o$B, D, K, 2, p, H)
v <- irf(f, id = id_chol(), horizon = H - 1, ci = "none")$estimate
rep_("T5 Cholesky orth IRF (exog=0 shocks, 13 horizons)", max(abs(ir$ort - v)) < 1e-9, "maxdiff", format(max(abs(ir$ort - v)), digits = 3))
cat("  impact matrix BEAR:\n"); print(round(ir$ort[, , 1], 4)); cat("  varm:\n"); print(round(v[, , 1], 4))

# reordering (BEAR Ordering: [ordered, rest]; P = chol(Sigma[ord,ord]); P[,backorder])
ord <- c(3, 1, 2); S2 <- o$Sigma[ord, ord]; Pn <- t(chol(S2)); bo <- order(ord)
# BEAR returns P(:,backorder) on a matrix whose rows are in the ORIGINAL ordering only after row back-ordering; we replicate the
# structural equivalent: D_orig = Pn[bo, bo-columns]
Dn <- matrix(0, K, K); Dn[ord, ord] <- Pn
ir2 <- bear_irfsim(o$B, Dn, K, 2, p, H)
v2 <- irf(f, id = id_chol(order = c("STN", "DOM_GDP", "DOM_CPI")), horizon = H - 1, ci = "none")$estimate
rep_("T6 Cholesky with reordering (STN first)", max(abs(ir2$ort - v2)) < 1e-9, "maxdiff", format(max(abs(ir2$ort - v2)), digits = 3))

# FEVD: BEAR fevd.m = cumulated squared orth IRF / total, gamma=1
fe <- array(0, c(K, K, H)); cs <- array(0, c(K, K, H))
for (h in 1:H) cs[, , h] <- (if (h == 1) 0 else cs[, , h - 1]) + ir$ort[, , h]^2
for (h in 1:H) fe[, , h] <- cs[, , h] / rowSums(cs[, , h])
fv <- fevd(f, id = id_chol(), horizon = H - 1, ci = "none")$estimate
rep_("T7 FEVD (13 horizons)", max(abs(fe - fv)) < 1e-9, "maxdiff", format(max(abs(fe - fv)), digits = 3))
rep_("T7 FEVD rows sum to 1", max(abs(apply(fv, c(1, 3), sum) - 1)) < 1e-10)

# HD: BEAR hdecomp.m: contribution_{jk,t} = sum_{l<=t} ortirf_{jk}(l) * shock_k(t-l+1), shocks = D^-1 u; baseline = Y - sum
T <- nrow(o$X); Sh <- t(solve(D) %*% t(o$E)); ir3 <- bear_irfsim(o$B, D, K, 2, p, T)
hdb <- array(0, c(T, K, K))
for (j in 1:K) for (k in 1:K) for (t in 1:T) hdb[t, j, k] <- sum(ir3$ort[j, k, 1:t] * rev(Sh[1:t, k]))
hv <- hd(f, id = id_chol())
rep_("T8 HD shock contributions", max(abs(hdb - hv$contributions)) < 1e-8, "maxdiff", format(max(abs(hdb - hv$contributions)), digits = 3))
base_b <- o$Y - apply(hdb, c(1, 2), sum)
rep_("T8 HD baseline (BEAR 'initial+deterministic' = Y - sum contribs)", max(abs(base_b - hv$baseline)) < 1e-8, "maxdiff", format(max(abs(base_b - hv$baseline)), digits = 3))
# exogenous/deterministic part separately (BEAR hdecomp_inc_exo): const + Oil effect is inside the baseline in both
cat("  (BEAR hdecomp_inc_exo splits baseline into deterministic + exogenous contributions; varm reports one baseline)\n")
