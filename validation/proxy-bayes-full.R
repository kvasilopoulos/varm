# Monte Carlo checks for varm_proxy_arw() and varm_proxy_gkr().
# Not run in tests. Usage: Rscript validation/proxy-bayes-full.R [reps] [T]
# DGP: K = 3, VAR(1), unit-variance shocks, impact matrix B0 (columns = shocks), proxy m = e3 + sd * v.
# Reports, over `reps` samples:
#   * ARW: coverage of the 68% / 90% posterior band for the impact response of y1, y2 to shock 3 (target y3 = 1)
#   * GKR: how often the true response to the set-identified shock 2 (sign: y1, y2 respond > 0 on impact)
#     lies in the set of posterior means / in the 90% robust credible region, and the mean width of each
#   * rate at which the single-prior (first accepted rotation) 90% band covers the truth
args <- as.numeric(commandArgs(TRUE))
reps <- if (length(args) >= 1) args[1] else 50
Tn <- if (length(args) >= 2) args[2] else 300
pkgload::load_all(".", quiet = TRUE)
B0 <- matrix(c(1, 0.3, 0.2, 0, 1, 0.4, 0.5, 0.2, 1), 3)
S <- matrix(0, 3, 3); S[1, 2] <- 1; S[2, 2] <- 1
tr3 <- B0[, 3] / B0[3, 3]
one <- function(seed) {
  set.seed(seed)
  e <- matrix(rnorm((Tn + 50) * 3), ncol = 3)
  y <- matrix(0, Tn + 50, 3)
  for (t in 2:(Tn + 50)) y[t, ] <- 0.4 * y[t - 1, ] + B0 %*% e[t, ]
  k <- 51:(Tn + 50)
  y <- y[k, ]; m <- e[k, 3] + rnorm(Tn, sd = 0.5)
  bf <- bvarm(y, p = 1, ndraw = 300, nburn = 100)
  a <- varm_proxy_arw(bf, m, target = 3, horizon = 1)
  g <- varm_proxy_gkr(bf, m, signs = S, shock = 2, horizon = 1, nq = 100, level = 0.9)
  c(
    arw68 = mean(a$bands[1, 1:2, 2] <= tr3[1:2] & tr3[1:2] <= a$bands[1, 1:2, 3]),
    arw90 = mean(a$bands[1, 1:2, 1] <= tr3[1:2] & tr3[1:2] <= a$bands[1, 1:2, 4]),
    gkr_set = mean(g$lower[1, 1:3] <= B0[1:3, 2] & B0[1:3, 2] <= g$upper[1, 1:3]),
    gkr_rcr = mean(g$region[1, 1:3, 1] <= B0[1:3, 2] & B0[1:3, 2] <= g$region[1, 1:3, 2]),
    gkr_bayes90 = mean(g$bands[1, 1:3, 1] <= B0[1:3, 2] & B0[1:3, 2] <= g$bands[1, 1:3, 2]),
    w_rcr = mean(g$region[1, 1:3, 2] - g$region[1, 1:3, 1]),
    w_bayes = mean(g$bands[1, 1:3, 2] - g$bands[1, 1:3, 1])
  )
}
res <- t(vapply(seq_len(reps), one, numeric(7)))
print(round(colMeans(res), 3))
