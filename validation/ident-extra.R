# Longer validation of proxy-SVAR, zero+sign and A-B restrictions. Run: Rscript validation/ident-extra.R
pkgload::load_all(".", quiet = TRUE)
set.seed(1)
Y <- as.matrix(yraw); fit <- varm(Y, p = 2)

## 1. proxy: simulate DGP with known shock, compare normalised IRF to truth; hand-coded gamma
t <- 2000; A <- matrix(c(.5, .1, 0, 0, .4, .1, .1, 0, .3), 3); P <- matrix(c(1, .3, .2, 0, 1, -.4, 0, 0, .8), 3)
e <- matrix(rnorm(3 * t), t); S <- matrix(0, t, 3); for (i in 2:t) S[i, ] <- A %*% S[i - 1, ] + P %*% e[i, ]
z <- e[, 2] + rnorm(t); f <- varm(S, 1)
r <- irf(f, id_proxy(z, f, target = 2), horizon = 4, ci = "none")$estimate[, 2, ]
truth <- sapply(0:4, function(h) { Ah <- diag(3); for (k in seq_len(h)) Ah <- Ah %*% A; (Ah %*% P[, 2]) / P[2, 2] })
cat("proxy max |est - truth| (T=2000):", max(abs(r - truth)), "\n")
for (ty in c("wild", "mbb")) {
  pd <- proxy_draws(f, z, 500, ty, target = 2); ir <- irf(pd$draws, pd$id, horizon = 4)
  inb <- ir$bands[, 2, , 1] <= truth & truth <= ir$bands[, 2, , 4]
  cat(ty, "90% band covers truth at", mean(inb), "of (variable,horizon) cells; F =", pd$fstat[["F"]], "\n")
}
## 2. coverage of the 90% proxy bootstrap band for the h=0 response of var 3, MC
cover <- function(ty) mean(replicate(60, {
  e <- matrix(rnorm(3 * 300), 300); S <- matrix(0, 300, 3); for (i in 2:300) S[i, ] <- A %*% S[i - 1, ] + P %*% e[i, ]
  f <- varm(S, 1); pd <- proxy_draws(f, e[, 2] + rnorm(300), 149, ty, target = 2)
  q <- quantile(irf(pd$draws, pd$id, horizon = 0)$draws[3, 2, 1, ], c(.05, .95))
  q[1] <= P[3, 2] / P[2, 2] && P[3, 2] / P[2, 2] <= q[2] }))
cat("MC coverage of nominal 90% band (h=0, var 3; 60 reps): mbb", cover("mbb"), " wild", cover("wild"), "\n")

## 3. zero+sign on real data: exactness and acceptance
d <- draws(fit, "bootstrap", 200)
Sg <- matrix(0, 3, 3); Sg[1, 1] <- -1; Sg[3, 1] <- 1; Sg[2, 1] <- 1   # a demand-like shock
Z <- matrix(0, 3, 3); Z[1, 2] <- 1
id <- id_zerosign(Sg, Z, horizon = 2, nrot = 5)
si <- varm:::structural_irf(d$B, d$Sigma, 3, 2, 12, id, colnames(Y))
cat("zero-sign accept rate:", si$accept_rate, " accepted:", length(si$draw), "\n")
cat("max zero violation:", max(abs(si$irf[1, 2, 1, ])), " max |PP'-Sigma|:",
    max(sapply(seq_along(si$draw), function(a) max(abs(si$P[, , a] %*% t(si$P[, , a]) - d$Sigma[, , si$draw[a]])))), "\n")
## 4. A-B vs vars::BQ
if (requireNamespace("vars", quietly = TRUE)) {
  v <- vars::VAR(Y, p = 2, type = "const"); bq <- vars::BQ(v)
  lt <- matrix(NA, 3, 3); lt[upper.tri(lt)] <- 0
  P1 <- irf(fit, id_restrictions(longrun = lt), horizon = 0, ci = "none")$estimate[, , 1]
  cat("longrun vs vars::BQ impact (abs, columns match up to sign):", max(abs(abs(P1) - abs(bq$B))), "\n")
}
