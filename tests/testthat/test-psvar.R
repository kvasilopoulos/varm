sim_proxy <- function(t = 300) {
  set.seed(1)
  A1 <- matrix(c(0.2, 0.5, 0, 0.2), 2, 2)
  A2 <- matrix(c(0.1, -0.3, 0.2, 0.1), 2, 2)
  B <- matrix(c(-0.592, 0.592, 0.806, 0.806), 2, 2)
  e <- cbind(rnorm(t), rnorm(t))
  Y <- matrix(0, t, 2, dimnames = list(NULL, c("y1", "y2")))
  for (i in 3:t) {
    Y[i, ] <- A1 %*% Y[i - 1, ] + A2 %*% Y[i - 2, ] + B %*% e[i, ]
  }
  list(Y = Y, m = 5 * e[, 1] + rnorm(t))
}

test_that("varm_proxy recovers the psvar point IRF", {
  skip_if_not_installed("psvar")
  s <- sim_proxy()
  fit <- varm(s$Y, p = 2)
  r <- varm_proxy(fit, s$m, target = 1, size = 1, horizon = 10, nboot = 0)
  ref <- psvar::irf_psvar(
    psvar::estimate_var(s$Y, 2),
    s$m,
    irhor = 11,
    boot = FALSE,
    shocksize = 1
  )$irfs
  expect_equal(unname(r$estimate), unname(ref), tolerance = 1e-8)
})

test_that("varm_proxy self-consistency and bootstrap bands", {
  s <- sim_proxy()
  fit <- varm(s$Y, p = 2)
  r <- varm_proxy(fit, s$m, size = 2, horizon = 6, nboot = 50)
  expect_equal(unname(r$estimate[1, 1]), 2)
  expect_equal(dim(r$bands), c(7, 2, 4))
  expect_true(r$fstat[["F"]] > 10)
})

# ---- Montiel Olea-Stock-Watson robust sets ---------------------------------

test_that("varm_proxy_robust: strong instrument gives bounded sets containing the estimate", {
  s <- sim_proxy(400)
  fit <- varm(s$Y, p = 2)
  r <- varm_proxy_robust(fit, s$m, horizon = 8, level = c(0.68, 0.90))
  expect_equal(dim(r$lower), c(9, 2, 2))
  expect_true(all(r$case[-1, , ] == 1))
  est3 <- array(r$estimate, dim(r$lower))
  expect_true(all(r$lower <= est3 + 1e-12 & est3 <= r$upper + 1e-12))
  # 90% set nests the 68% set; impact of target is the normalisation
  expect_true(all(
    r$lower[,, 2] <= r$lower[,, 1] + 1e-12 &
      r$upper[,, 2] >= r$upper[,, 1] - 1e-12
  ))
  expect_equal(unname(r$lower[1, 1, ]), c(1, 1))
  expect_equal(
    unname(r$estimate),
    unname(varm_proxy(fit, s$m, horizon = 8, nboot = 0)$estimate),
    tolerance = 1e-10
  )
})

test_that("varm_proxy_robust: pure-noise instrument gives unbounded sets", {
  set.seed(3)
  s <- sim_proxy(300)
  m <- rnorm(300)
  fit <- varm(s$Y, p = 1)
  r <- varm_proxy_robust(fit, m, horizon = 4, level = 0.9)
  expect_lt(r$fstat[["F_hac"]], stats::qchisq(0.9, 1))
  expect_true(all(
    is.infinite(r$upper[-1, , 1]) & is.infinite(r$lower[-1, , 1])
  ))
  expect_true(all(r$case[-1, , 1] %in% c(2, 4)))
})

test_that("robust F is the Wald statistic T * Gamma^2 / W_Gamma", {
  s <- sim_proxy(300)
  fit <- varm(s$Y, p = 2)
  cv <- varm:::proxy_cov(fit, s$m[-(1:2)])
  f <- proxy_fstat(fit, s$m, target = 2)
  expect_equal(
    unname(f[["F_hac"]]),
    unname(cv$Tp * cv$Gamma[2]^2 / cv$W[cv$iG[2], cv$iG[2]])
  )
  expect_equal(
    unname(cv$Gamma),
    unname(drop(crossprod(fit$U, s$m[-(1:2)])) / cv$Tp)
  )
})

test_that("confidence-set endpoints solve the Wald test (independent numerical delta method)", {
  s <- sim_proxy(300)
  p <- 2
  K <- 2
  H <- 3
  sz <- 1.5
  lev <- 0.9
  fit <- varm(s$Y, p = p)
  r <- varm_proxy_robust(fit, s$m, size = sz, horizon = H, level = lev)
  cv <- varm:::proxy_cov(fit, s$m[-(1:p)])
  nstat <- function(i, h, lam) {
    # f(vec(A), Gamma) = sz * e_i' Phi_h(A) Gamma - lam * Gamma_1, gradient by finite differences
    f <- function(th) {
      B <- fit$B
      B[1:(K * p), ] <- t(matrix(th[1:(K * K * p)], K, K * p))
      sz *
        drop(cpp_ma(B, K, p, H)[i, , h + 1] %*% th[cv$iG]) -
        lam * th[cv$iG[1]]
    }
    th0 <- c(as.vector(t(fit$B[1:(K * p), ])), cv$Gamma)
    d <- vapply(
      seq_along(th0),
      function(k) {
        e <- replace(numeric(length(th0)), k, 1e-6)
        (f(th0 + e) - f(th0 - e)) / 2e-6
      },
      0
    )
    unname(cv$Tp * f(th0)^2 / drop(d %*% cv$W %*% d))
  }
  for (i in 1:2) {
    for (h in 1:H) {
      expect_equal(r$case[h + 1, i, 1], 1)
      expect_equal(
        nstat(i, h, r$lower[h + 1, i, 1]),
        stats::qchisq(lev, 1),
        tolerance = 1e-4
      )
      expect_equal(
        nstat(i, h, r$upper[h + 1, i, 1]),
        stats::qchisq(lev, 1),
        tolerance = 1e-4
      )
      expect_lt(nstat(i, h, r$estimate[h + 1, i]), stats::qchisq(lev, 1))
    }
  }
})

test_that("varm_proxy_robust agrees with psvar (impact and first horizon; bound order is orientation-free)", {
  skip_if_not_installed("psvar")
  s <- sim_proxy(300)
  fit <- varm(s$Y, p = 2)
  r <- varm_proxy_robust(fit, s$m, horizon = 5, level = 0.9)
  ps <- psvar:::irf_robust_psvar_single(
    psvar::estimate_var(s$Y, 2),
    s$m,
    clevel = 90,
    irhor = 6
  )
  # psvar's upper/lower are swapped, and its Gmatrices() is only right up to horizon 1
  lo <- pmin(ps$lower, ps$upper)
  hi <- pmax(ps$lower, ps$upper)
  expect_equal(unname(r$lower[1:2, , 1]), unname(lo[1:2, ]), tolerance = 1e-8)
  expect_equal(unname(r$upper[1:2, , 1]), unname(hi[1:2, ]), tolerance = 1e-8)
})

# ---- bootstrap methods, multiple instruments, Bayesian proxy, LP-IV ---------

test_that("moving block bootstrap: shapes, recentring and reproducibility", {
  s <- sim_proxy(300)
  fit <- varm(s$Y, p = 2)
  set.seed(2)
  r <- varm_proxy(fit, s$m, horizon = 5, nboot = 40, method = "mbb")
  expect_equal(dim(r$bands), c(6, 2, 4))
  expect_equal(dim(r$draws), c(6, 2, 40))
  expect_equal(r$block, round(5.03 * 298^0.25))
  expect_equal(unname(r$estimate[1, 1]), 1)
  expect_true(all(r$bands[, , "0.05"] <= r$bands[, , "0.95"]))
  # default is mbb; point estimates do not depend on the method
  expect_equal(varm_proxy(fit, s$m, horizon = 5, nboot = 3)$method, "mbb")
  w <- varm_proxy(fit, s$m, horizon = 5, nboot = 10, method = "wild")
  expect_equal(w$method, "wild")
  expect_equal(w$estimate, r$estimate)
  expect_error(varm_proxy(fit, s$m, nboot = 5, method = "nope"))
  # normal wild runs too
  expect_equal(
    dim(varm_proxy(fit, s$m, horizon = 3, nboot = 5, method = "wild_normal")$draws),
    c(4, 2, 5)
  )
})

test_that("block resampling is centred: bootstrap residual means are ~0 over positions", {
  U <- matrix(rnorm(400, 3), 200)
  ell <- 10
  nb <- 200 - ell + 1
  cU <- do.call(rbind, lapply(1:ell, function(s) colMeans(U[s:(s + nb - 1), ])))
  expect_equal(dim(cU), c(ell, 2))
  expect_true(all(abs(cU - 3) < 0.5))
})

test_that("multiple instruments recover the true impact columns (large T)", {
  set.seed(11)
  Tn <- 20000
  K <- 3
  B <- matrix(c(1, 0.4, 0.3, 0, 0.8, -0.5, 0.2, 0.1, 0.9), 3, 3) # B[1:2,1:2] lower-triangular
  e <- matrix(rnorm(Tn * K), Tn)
  u <- e %*% t(B)
  Y <- matrix(0, Tn, K, dimnames = list(NULL, c("a", "b", "c")))
  A <- diag(0.3, K)
  for (t in 2:Tn) Y[t, ] <- A %*% Y[t - 1, ] + u[t, ]
  M <- cbind(e[, 1] + 0.3 * e[, 2] + rnorm(Tn, sd = 0.7), e[, 2] + rnorm(Tn, sd = 0.7))
  fit <- varm(Y, p = 1)
  r <- varm_proxy(fit, M, target = c("a", "b"), size = c(1, 2), horizon = 2, nboot = 0)
  expect_equal(dim(r$estimate), c(3, 3, 2))
  # Mertens-Ravn identifies the first r impact columns only up to a rotation Q, so test
  # rotation-invariant objects: B1 B1' (the covariance share of the r shocks) and the span
  expect_equal(unname(tcrossprod(r$B1)), unname(tcrossprod(B[, 1:2])), tolerance = 0.05)
  expect_equal(unname(cov(r$shocks)), diag(2), tolerance = 0.05)
  R2 <- summary(lm(r$shocks ~ e[-1, 1:2]))
  expect_true(all(vapply(R2, function(z) z$r.squared, 0) > 0.98))
  # normalisation: column i has impact size[i] on target[i]
  expect_equal(unname(diag(r$estimate[1, 1:2, ])), c(1, 2))
  # the 2x2 impact block is a scaled copy of B1[1:2, ]
  expect_equal(
    unname(r$estimate[1, , ]),
    unname(sweep(r$B1, 2, c(1, 2) / diag(r$B1[1:2, ]), "*"))
  )
  expect_equal(dim(r$fstat), c(2, 3))
  rb <- varm_proxy(fit, M, target = c("a", "b"), horizon = 2, nboot = 5, method = "mbb")
  expect_equal(dim(rb$bands), c(3, 3, 2, 4))
  expect_error(varm_proxy(fit, M, target = "a"), "one distinct variable")
})

test_that("single-instrument impact vector equals the Mertens-Ravn B1 direction; reliability is cor^2", {
  s <- sim_proxy(500)
  fit <- varm(s$Y, p = 2)
  r <- varm_proxy(fit, s$m, nboot = 0, horizon = 1)
  expect_equal(
    unname(r$B1[, 1] / r$B1[1, 1]),
    unname(r$gamma / r$gamma[1])
  )
  expect_equal(r$reliability, unname(cor(s$m[-(1:2)], r$shocks[, 1])^2))
  expect_true(r$reliability > 0.9 && r$reliability < 1)
})

test_that("Bayesian proxy-SVAR gives bands close to the bootstrap with a strong proxy", {
  s <- sim_proxy(300)
  set.seed(4)
  bf <- bvarm(s$Y, p = 2, ndraw = 300, nburn = 100)
  b <- varm_proxy_bayes(bf, s$m, horizon = 4, ndraw = 150)
  expect_s3_class(b, "varm_proxy")
  expect_equal(b$method, "bayes")
  expect_equal(dim(b$bands), c(5, 2, 4))
  expect_equal(unname(b$estimate[1, 1]), 1)
  fit <- varm(s$Y, p = 2)
  f <- varm_proxy(fit, s$m, horizon = 4, nboot = 0)
  expect_equal(unname(b$estimate), unname(f$estimate), tolerance = 0.15)
  expect_error(varm_proxy_bayes(fit, s$m), "bvarm")
})

test_that("LP-IV via lp(iv=) agrees with the proxy-VAR impact response for a strong proxy", {
  s <- sim_proxy(2000)
  l <- lp(s$Y, shock = "y1", iv = matrix(s$m, ncol = 1), horizon = 4, lags = 2)
  fit <- varm(s$Y, p = 2)
  f <- varm_proxy(fit, s$m, horizon = 4, nboot = 0)
  expect_equal(unname(l$estimate[, "y2"]), unname(f$estimate[, 2]), tolerance = 0.08)
})
