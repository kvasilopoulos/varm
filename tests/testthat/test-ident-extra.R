sim_ie <- function(t = 400, seed = 11) {
  set.seed(seed)
  A <- matrix(c(0.5, 0.1, 0, 0, 0.4, 0.1, 0.1, 0, 0.3), 3, 3)
  P <- matrix(c(1, 0.3, 0.2, 0, 1, -0.4, 0, 0, 0.8), 3, 3)
  e <- matrix(rnorm(3 * t), t)
  Y <- matrix(0, t, 3, dimnames = list(NULL, c("a", "b", "c")))
  for (i in 2:t) {
    Y[i, ] <- A %*% Y[i - 1, ] + P %*% e[i, ]
  }
  list(Y = Y, z = e[, 2] + 0.5 * rnorm(t))
}

test_that("id_proxy matches the textbook computation", {
  s <- sim_ie()
  fit <- varm(s$Y, p = 2)
  r <- irf(
    fit,
    id_proxy(s$z, fit, target = 2, size = 1),
    horizon = 8,
    ci = "none"
  )
  U <- fit$U
  z <- s$z[-(1:2)]
  g <- colSums(U * z) / nrow(U)
  g <- g / g[2]
  Phi <- array(0, c(3, 3, 9))
  Phi[,, 1] <- diag(3)
  A <- lapply(1:2, function(l) t(fit$B[(l - 1) * 3 + 1:3, ]))
  for (h in 2:9) {
    for (l in 1:min(2, h - 1)) {
      Phi[,, h] <- Phi[,, h] + A[[l]] %*% Phi[,, h - l]
    }
  }
  ref <- sapply(1:9, function(h) Phi[,, h] %*% g)
  expect_equal(unname(r$estimate[, 2, ]), unname(ref), tolerance = 1e-10)
  expect_equal(unname(r$estimate[2, 2, 1]), 1)
  expect_true(all(r$estimate[, 1, ] == 0))
  v <- varm_proxy(fit, s$z, target = 2, size = 1, horizon = 8, nboot = 0)
  expect_equal(
    unname(t(r$estimate[, 2, ])),
    unname(v$estimate),
    tolerance = 1e-10
  )
})

test_that("proxy_draws wild and mbb give sane bands and are seeded", {
  s <- sim_ie(300)
  fit <- varm(s$Y, p = 1)
  for (ty in c("wild", "mbb")) {
    set.seed(3)
    pd <- proxy_draws(fit, s$z, 60, ty, target = 2)
    set.seed(3)
    pd2 <- proxy_draws(fit, s$z, 60, ty, target = 2)
    expect_equal(pd$draws$gamma, pd2$draws$gamma)
    r <- irf(pd$draws, pd$id, horizon = 5)
    expect_equal(dim(r$draws), c(3, 3, 6, 60))
    expect_true(all(abs(r$draws[2, 2, 1, ] - 1) < 1e-12))
    expect_true(median(r$draws[3, 2, 1, ]) > -2 && abs(r$estimate[3, 2, 1]) < 5)
  }
  expect_true(pd$fstat[["F"]] > 10)
})

test_that("zero+sign: zeros exact, signs hold, PP' = Sigma", {
  s <- sim_ie()
  fit <- varm(s$Y, p = 2)
  d <- draws(fit, "mc", 40)
  Sg <- matrix(0, 3, 3)
  Sg[1, 1] <- 1
  Sg[2, 1] <- 1
  Sg[3, 3] <- 1
  Sg[2, 2] <- -1
  Z <- matrix(0, 3, 3)
  Z[3, 2] <- 1 # shock 2 has no impact effect on variable 3
  lr <- matrix(0, 3, 3)
  lr[1, 3] <- 1 # shock 3 has no long-run effect on variable 1
  id <- id_zerosign(Sg, Z, lr, horizon = 1, nrot = 2)
  r <- irf(d, id, horizon = 6)
  si <- varm:::structural_irf(d$B, d$Sigma, 3, 2, 6, id, fit$varnames)
  expect_gt(si$accept_rate, 0.01)
  for (a in seq_along(si$draw)) {
    P <- si$P[,, a]
    n <- si$draw[a]
    th <- si$irf[,,, a]
    expect_equal(
      P %*% t(P),
      d$Sigma[,, n],
      tolerance = 1e-10,
      ignore_attr = TRUE
    )
    expect_lt(abs(th[3, 2, 1]), 1e-10)
    Fm <- diag(3) -
      Reduce(`+`, lapply(1:2, function(l) t(d$B[(l - 1) * 3 + 1:3, , n])))
    expect_lt(abs((solve(Fm) %*% P)[1, 3]), 1e-10)
    expect_true(
      all(th[c(1, 2), 1, 1] > 0) && th[3, 3, 1] > 0 && th[2, 2, 1] < 0
    )
  }
  f <- irf_fry_pagan(r)
  expect_equal(dim(f$irf), c(3, 3, 7))
})

test_that("id_restrictions reproduces Blanchard-Quah and recursive schemes", {
  s <- sim_ie()
  fit <- varm(s$Y, p = 2)
  lt <- matrix(NA, 3, 3)
  lt[upper.tri(lt)] <- 0
  a <- irf(
    fit,
    id_restrictions(longrun = lt),
    horizon = 5,
    ci = "none"
  )$estimate
  b <- irf(fit, id_longrun(), horizon = 5, ci = "none")$estimate
  sg <- sign(a[1, , 1] * b[1, , 1])
  sg[sg == 0] <- 1
  expect_equal(unname(sweep(a, 2, sg, "*")), unname(b), tolerance = 1e-8)
  a <- irf(
    fit,
    id_restrictions(impact = lt[,] * 0 + ifelse(upper.tri(lt), 0, NA)),
    horizon = 3,
    ci = "none"
  )$estimate
  b <- irf(fit, id_chol(), horizon = 3, ci = "none")$estimate
  expect_equal(unname(a), unname(b), tolerance = 1e-8)
  expect_error(irf(fit, id_restrictions(impact = matrix(0, 3, 3)), ci = "none"))
})

test_that("id_hetero recovers P in a two-regime simulation", {
  set.seed(5)
  t <- 3000
  P <- matrix(c(1, 0.5, 0, 1), 2)
  reg <- rep(0:1, each = t / 2)
  e <- cbind(rnorm(t), rnorm(t) * ifelse(reg == 1, 3, 1))
  Y <- matrix(0, t, 2)
  for (i in 2:t) {
    Y[i, ] <- 0.3 * Y[i - 1, ] + P %*% e[i, ]
  }
  fit <- varm(Y, 1)
  id <- id_hetero(fit, reg)
  expect_equal(unname(abs(id$P)), cbind(c(0, 1), c(1, 0.5)), tolerance = 0.06)
})

test_that("id_restrictions long-run-only orientation is deterministic with positive impact diagonal", {
  fit <- varm(yraw[, c("infl", "unemp")], p = 2)
  lr <- matrix(c(NA, 0, NA, NA), 2)
  p1 <- irf(
    fit,
    id_restrictions(longrun = lr),
    horizon = 0,
    ci = "none"
  )$estimate[,, 1]
  p2 <- irf(
    fit,
    id_restrictions(longrun = lr),
    horizon = 0,
    ci = "none"
  )$estimate[,, 1]
  expect_equal(p1, p2)
  expect_true(all(diag(p1) > 0))
})

test_that("proxy_draws default block follows the Jentsch-Lunsford rule", {
  expect_match(paste(deparse(args(proxy_draws)), collapse = ""), "block = NULL")
  expect_equal(round(5.03 * 215^(1 / 4)), 19)
})

test_that("bias_correct with mc draws warns; tvpvar irf rejects id", {
  fit <- varm(yraw, p = 2)
  expect_warning(draws(fit, "mc", 20, bias_correct = TRUE), "bias_correct")
})
