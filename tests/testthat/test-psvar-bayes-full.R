# Full-posterior Bayesian proxy-SVAR (ARW-type) and GKR robust-Bayes bounds

sim_pb <- function(Tn = 500, noise = 0.3, seed = 11) {
  set.seed(seed)
  B0 <- matrix(c(1, 0.3, 0.2, 0, 1, 0.4, 0.5, 0.2, 1), 3) # columns = shocks
  e <- matrix(rnorm(Tn * 3), Tn)
  y <- matrix(0, Tn, 3, dimnames = list(NULL, c("a", "b", "c")))
  for (t in 2:Tn) {
    y[t, ] <- 0.4 * y[t - 1, ] + B0 %*% e[t, ]
  }
  list(
    y = y,
    e = e,
    B0 = B0,
    m3 = e[, 3] + rnorm(Tn, sd = noise),
    m23 = cbind(e[, 2] + rnorm(Tn, sd = noise), e[, 3] + rnorm(Tn, sd = noise))
  )
}
s <- sim_pb()
bf <- bvarm(s$y, p = 1, ndraw = 300, nburn = 100)

test_that("rotation sampler satisfies exogeneity and Sigma = P P'", {
  fs <- proxy_first_stage(bf, s$m3)
  fd <- proxy_first_stage_draw(fs, 5)
  o <- cpp_proxy_rot(
    bf$draws$B[,, 1:5],
    bf$draws$Sigma[,, 1:5],
    fd$D,
    1L,
    0L,
    array(0, c(3, 3, 1)),
    3L,
    3L,
    100L,
    0L,
    -1L,
    0
  )
  expect_length(o$draw, 15L)
  for (a in 1:15) {
    P <- matrix(o$irf[,, a], 3, 3)
    n <- o$draw[a]
    expect_equal(
      P %*% t(P),
      matrix(bf$draws$Sigma[,, n], 3, 3),
      tolerance = 1e-10
    )
    expect_lt(max(abs(matrix(fd$D[,, n], 1, 3) %*% P[, 1:2])), 1e-9)
  }
})

test_that("varm_proxy_arw: strong proxy recovers the true impact column and matches varm_proxy", {
  a <- varm_proxy_arw(bf, s$m3, target = 3, horizon = 4)
  true <- s$B0[, 3] / s$B0[3, 3]
  expect_equal(unname(a$estimate[1, ]), true, tolerance = 0.08)
  expect_true(all(
    a$bands[1, , 1] <= true + 1e-8 & true <= a$bands[1, , 4] + 1e-8
  ))
  v <- varm_proxy(varm(s$y, p = 1), s$m3, target = 3, horizon = 4, nboot = 0)
  expect_equal(unname(a$estimate), unname(v$estimate), tolerance = 0.05)
  # very informative instrument: posterior collapses on the frequentist direction
  s2 <- sim_pb(noise = 0.02)
  b2 <- bvarm(s2$y, p = 1, ndraw = 200, nburn = 100)
  a2 <- varm_proxy_arw(b2, s2$m3, target = 3, horizon = 4)
  v2 <- varm_proxy(varm(s2$y, p = 1), s2$m3, target = 3, horizon = 4, nboot = 0)
  expect_equal(
    unname(a2$estimate[1, ]),
    unname(v2$estimate[1, ]),
    tolerance = 0.02
  )
  expect_equal(unname(a2$estimate), unname(v2$estimate), tolerance = 0.06) # dynamics differ by the Minnesota prior
  expect_s3_class(a, "varm_proxy")
})

test_that("id_proxy_bayes plugs into irf() and reliability truncation filters draws", {
  ir <- irf(bf, id_proxy_bayes(bf, s$m3, target = 3), horizon = 4)
  expect_s3_class(ir, "varm_dyn")
  expect_equal(unname(ir$estimate[1, 3, 1]), 0.5, tolerance = 0.1) # impact of shock 3 on a
  expect_true(all(ir$estimate[, 1:2, ] == 0)) # non-proxy shocks are not identified
  a_all <- varm_proxy_arw(bf, s$m3, target = 3, horizon = 2, reliability = 0)
  a_hi <- varm_proxy_arw(bf, s$m3, target = 3, horizon = 2, reliability = 0.92)
  expect_lt(a_hi$nacc, a_all$nacc)
  expect_error(
    varm_proxy_arw(bf, s$m3, target = 3, horizon = 2, reliability = 0.999),
    "reliability"
  )
})

test_that("two proxies for two shocks: sign restrictions are respected by every draw", {
  sg <- matrix(0, 3, 2)
  sg[2, 1] <- 1 # shock 2 raises b
  sg[3, 2] <- 1 # shock 3 raises c
  a <- varm_proxy_arw(
    bf,
    s$m23,
    target = c(2, 3),
    size = NULL,
    signs = sg,
    horizon = 2,
    nrot = 2
  )
  expect_equal(dim(a$draws)[3], 2L)
  expect_true(all(a$draws[1, 2, 1, ] > 0))
  expect_true(all(a$draws[1, 3, 2, ] > 0))
  expect_lt(a$accept_rate, 1)
  # with r = 2 the signs only restrict the rotation: the true columns are inside the posterior range
  for (j in 1:2) {
    rg <- apply(a$draws[1, , j, ], 1, range)
    expect_true(all(
      rg[1, ] - 0.05 <= s$B0[, j + 1] & s$B0[, j + 1] <= rg[2, ] + 0.05
    ))
  }
})

test_that("GKR: bounds contain the truth, mean of single-prior draws lies in the set, point-ID shock collapses", {
  S <- matrix(0, 3, 3)
  S[1, 2] <- 1 # shock 2 (non-proxy) raises a and b on impact
  S[2, 2] <- 1
  g <- varm_proxy_gkr(bf, s$m3, signs = S, shock = 2, horizon = 3, nq = 150)
  expect_s3_class(g, "varm_proxy_gkr")
  expect_equal(g$prob_nonempty, 1)
  expect_true(all(g$lower <= g$upper + 1e-12))
  # identified set at the truth contains the true (sd-scaled) impact response of shock 2
  expect_true(all(
    g$lower[1, ] - 0.1 <= s$B0[, 2] & s$B0[, 2] <= g$upper[1, ] + 0.1
  ))
  # set of posterior means contains the posterior mean of any single prior (here: first accepted rotation)
  mdraw <- apply(g$draws, 1:2, mean)
  expect_true(all(g$lower - 1e-10 <= mdraw & mdraw <= g$upper + 1e-10))
  # robust credible regions widen with the level (columns: lower90, lower68, upper68, upper90)
  expect_true(all(
    g$region[,, 4] - g$region[,, 1] >= g$region[,, 3] - g$region[,, 2] - 1e-9
  ))
  # unrestricted shock (K) is point identified: bounds coincide and equal the ARW / varm_proxy answer
  g3 <- varm_proxy_gkr(bf, s$m3, shock = 3, target = 3, size = 1, horizon = 3)
  expect_equal(g3$lower, g3$upper, tolerance = 1e-8)
  v <- varm_proxy(varm(s$y, p = 1), s$m3, target = 3, horizon = 3, nboot = 0)
  expect_equal(unname(g3$lower), unname(v$estimate), tolerance = 0.05)
})

test_that("GKR with r = 2 proxies: first shock is point identified; incompatible signs give an error", {
  g <- varm_proxy_gkr(bf, s$m23, shock = 1, horizon = 2, nq = 20)
  expect_equal(g$lower, g$upper, tolerance = 1e-8)
  # the proxy shock is point identified up to sign (set by diag(A0) >= 0): exactly one of +c / -c is non-empty
  empty <- vapply(
    c(-1, 1),
    function(sgn) {
      Sx <- matrix(0, 3, 3)
      Sx[3, 3] <- sgn
      inherits(
        try(
          varm_proxy_gkr(
            bf,
            s$m3,
            signs = Sx,
            shock = 3,
            horizon = 2,
            nq = 20,
            maxtry = 50
          ),
          silent = TRUE
        ),
        "try-error"
      )
    },
    NA
  )
  expect_equal(sum(empty), 1L)
})
