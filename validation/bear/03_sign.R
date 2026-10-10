# BEAR irfres.m sign-restriction sampler (sequential column Q, checksignres.m) vs varm id_sign, at FIXED OLS (B, Sigma).
# Both are Monte-Carlo samplers of the same target (Haar rotations conditional on the sign restrictions), so we compare
# IRF quantiles across accepted rotations (MC error ~ 1/sqrt(N)).
source("validation/bear/common.R")
set.seed(7)
p <- 4; d <- load_bear_data(p); Y <- d$Y; X <- d$X; K <- 3; H <- 5; Hs <- 4
f <- varm(Y, p, const = TRUE, exog = X); o <- bear_olsvar(Y, X, TRUE, p)
Dl <- t(chol(o$Sigma)); ir <- bear_irfsim(o$B, Dl, K, 2, p, H)$ort   # [resp, shock, h], h=1..H
# restrictions (resp x shock): shock1 demand (+,+,+) shock2 supply (GDP+,CPI-) shock3 policy (GDP-,CPI-,STN+)
S <- matrix(0, K, K); S[, 1] <- c(1, 1, 1); S[, 2] <- c(1, -1, 0); S[, 3] <- c(-1, -1, 1)
stk <- lapply(1:K, function(j) NULL)
chk1 <- function(j, f) { s <- as.vector(S[, j]); idx <- s != 0
  if (!any(idx)) return(1); v <- rep(s[idx], Hs) * unlist(lapply(1:Hs, function(h) f[idx, h]))  # f: K x Hs response of shock candidate
  if (all(v >= 0)) 1 else if (all(-v >= 0)) -1 else 0 }
bear_sign <- function(n) {
  out <- array(0, c(K, K, H, n)); tries <- 0; got <- 0
  while (got < n) {
    tries <- tries + 1; Q <- matrix(0, K, 0); ok <- TRUE
    for (j in 1:K) {
      x <- rnorm(K)
      q <- if (ncol(Q) == 0) x / sqrt(sum(x^2)) else { N <- MASS::Null(Q); y <- drop(N %*% (t(N) %*% x)); y / sqrt(sum(y^2)) }
      fj <- sapply(1:Hs, function(h) ir[, , h] %*% q)  # K x Hs
      s <- chk1(j, fj); if (s == 0) { ok <- FALSE; break }; Q <- cbind(Q, s * q)
    }
    if (ok) { got <- got + 1; for (h in 1:H) out[, , h, got] <- ir[, , h] %*% Q }
  }
  list(irf = out, rate = got / tries)
}
n <- 3000; b <- bear_sign(n)
id <- id_sign(S, horizon = Hs, nrot = n, maxtry = 1e6)
v <- irf(f, id = id, horizon = H - 1, ci = "none", level = 0.68)
a4 <- v$draws; cat("varm accepted:", dim(a4)[4], " accept rate:", format(v$accept_rate, digits = 3), " BEAR-replica rate:", format(b$rate, digits = 3), "\n")
qs <- function(a, h, i, j) quantile(a[i, j, h, ], c(.16, .5, .84))
mx <- 0; mxz <- 0
for (h in c(1, 3, 5)) for (i in 1:K) for (j in 1:K) {
  qb <- qs(b$irf, h, i, j); qv <- qs(a4, h, i, j); sc <- sd(b$irf[i, j, h, ])
  mx <- max(mx, max(abs(qb - qv)) / sc)
}
mb <- apply(b$irf, 1:3, mean); mv <- apply(a4, 1:3, mean); sdv <- apply(b$irf, 1:3, sd)
z <- max(abs(mb - mv)[, , c(1, 3, 5)] / (sdv[, , c(1, 3, 5)] * sqrt(2 / n)))
cat("max quantile diff / sd:", round(mx, 3), "  max |z| of mean diff:", round(z, 2), "\n")
rep_("T9 sign-restriction IRF distribution (3 shocks, 4 periods) varm vs BEAR-replica", z < 4.5 && mx < 0.2)
viol <- 0; for (r in seq_len(dim(a4)[4])) for (j in 1:K) for (h in 1:Hs) viol <- viol + sum(S[, j] * a4[, j, h, r] < 0)
rep_("T9b restrictions satisfied by every varm draw", viol == 0, "violations:", viol)

## ---- T10: BEAR-shipped sign-restriction output (results_ols.mat D_record, 1001 draws) vs varm ----
# BEAR OLS run: sign table Demand(+,+,+) Supply(GDP+,CPI-,.) Money(GDP+,CPI+,STN-) at impact + relative-magnitude restriction
# (GDP response to Demand >= GDP response to Supply), parameters drawn by olsirft4.m: Sigma^-1 ~ W(Sigma_hat^-1/T, T), beta ~ N(betahat, Sigma (x) (X'X)^-1).
set.seed(11)
o <- bear_mat("results_ols.mat", c("data", "D_record")); Yb <- o$data[, 1:3]; fo <- varm(Yb, 4)
S2 <- matrix(0, 3, 3); S2[, 1] <- c(1, 1, 1); S2[, 2] <- c(1, -1, 0); S2[, 3] <- c(1, 1, -1)
res <- irf(fo, id = id_sign(S2, horizon = 1, nrot = 1, maxtry = 500), horizon = 0, ci = "mc", nboot = 60000, level = .9)
a <- res$draws[, , 1, ]                      # [resp, shock, N]
keep <- a[1, 1, ] >= a[1, 2, ]
cat("varm: accepted rotations:", dim(a)[3], " after relative-magnitude filter:", sum(keep), "\n")
vb <- t(matrix(a[, , keep], 9)); bb <- t(o$D_record)   # BEAR D_record column-major: [resp, shock]
qn <- c(.05, .5, .95)
tab <- sapply(1:9, function(i) c(bear = quantile(bb[, i], qn), varm = quantile(vb[, i], qn)))
rownames(tab) <- c(paste0("BEAR", qn), paste0("varm", qn)); colnames(tab) <- paste0(rep(c("GDP", "CPI", "STN"), 3), "<-", rep(c("Dem", "Sup", "Mon"), each = 3))
print(round(tab, 3))
mad <- max(abs(tab[4:6, ] - tab[1:3, ])); scl <- apply(bb, 2, sd)
zmed <- max(abs(tab["varm0.5", ] - tab["BEAR0.5", ]) / scl)
cat("max |median diff| / sd(BEAR draws):", round(zmed, 3), "   max quantile diff:", round(mad, 3), "\n")
rep_("T10 BEAR stored sign-restricted D_record vs varm (mc draws + sign + relmagn filter), medians within 0.15 sd", zmed < 0.15)

## ---- T13 BEAR find_medianmodel.m (standardised median-target, 'MM=1') vs varm Fry-Pagan (irf_fry_pagan / varm:::fry_pagan) on the same draws ----
bear_mm <- function(a4) { d <- dim(a4); N <- d[4]; tot <- numeric(N)
  for (i in 1:d[1]) for (j in 1:d[2]) for (h in 1:d[3]) { x <- a4[i, j, h, ]; tot <- tot + ((x - median(x)) / sd(x))^2 }
  which.min(tot) }
ib <- bear_mm(a4); iv <- varm:::fry_pagan(a4)
rep_("T13 median-model (Fry-Pagan) index, BEAR find_medianmodel vs varm fry_pagan", ib == iv, "BEAR", ib, "varm", iv)
