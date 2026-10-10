# Statistical validation of the Bayesian engines (fast versions; full study: validation/bayes.R)

# K = 1 triangular posterior against 1-D quadrature of the exact semi-conjugate posterior
test_that("triangular homoskedastic sampler (K = 1) matches the exact posterior", {
  set.seed(1)
  y <- matrix(as.numeric(yraw[[1]]), ncol = 1)
  p <- 2
  pr <- prior_triangular("minnesota", lambda1 = 0.3)
  f <- bvarm(y, p = p, prior = pr, ndraw = 20000, nburn = 300)
  X <- cpp_var_ols(y, p, make_det(nrow(y), TRUE, FALSE, NULL))$X
  yy <- y[-(1:p), 1]
  sc <- sqrt(ar_resid_var(y, p))
  v <- c(sc^2 / sc^2 / (1:p)^pr$alpha * pr$lambda1^2, (pr$lambda3 * sc)^2)
  g0 <- c(1, rep(0, p - 1), 0)
  r <- yy - X %*% g0
  ev <- eigen(X %*% (v * t(X)), symmetric = TRUE)
  dk <- pmax(ev$values, 0)
  ur2 <- crossprod(ev$vectors, r)[, 1]^2
  ls <- seq(log(sc^2) - 6, log(sc^2) + 4, length.out = 4000)
  s2 <- exp(ls)
  lw <- vapply(
    s2,
    function(s) -0.5 * sum(log(s + dk)) - 0.5 * sum(ur2 / (s + dk)),
    0
  ) -
    (pr$hom[1] + 1) * ls -
    pr$hom[2] * sc^2 / s2 +
    ls
  w <- exp(lw - max(lw))
  w <- w / sum(w)
  mom <- vapply(
    seq_along(s2),
    function(k) {
      P <- crossprod(X) / s2[k] + diag(1 / v)
      m <- solve(P, crossprod(X, yy) / s2[k] + g0 / v)
      c(m, diag(solve(P)))
    },
    numeric(2 * ncol(X))
  )
  mu <- mom[seq_len(ncol(X)), ] %*% w
  ES2 <- sum(w * s2)
  ev_ <- sum(w * s2)
  vs2 <- sum(w * s2^2) - ev_^2
  expect_lt(
    max(
      abs(rowMeans(matrix(f$draws$B, ncol(X))) - mu) /
        (apply(matrix(f$draws$B, ncol(X)), 1, sd) / sqrt(20000))
    ),
    5
  )
  expect_equal(mean(f$draws$Sigma), ES2, tolerance = 0.02)
  expect_equal(var(as.vector(f$draws$Sigma)), vs2, tolerance = 0.08)
})

test_that("loose triangular prior reproduces OLS", {
  set.seed(2)
  Y <- as.matrix(yraw[, 1:3])
  f <- bvarm(
    Y,
    p = 2,
    prior = prior_triangular(
      "minnesota",
      lambda1 = 1e3,
      lambda2 = 1,
      lambda3 = 1e3,
      lambda5 = 1e3,
      hom = c(0.01, 0.01)
    ),
    ndraw = 4000,
    nburn = 200
  )
  expect_lt(max(abs(coef(f) - varm(Y, p = 2)$B)), 0.03)
})

test_that("KSC / FFBS volatility update passes the Geweke joint-distribution test", {
  set.seed(3)
  g <- cpp_sv_geweke(T = 10, nsim = 150000, a_h = 5, b_h = 0.2, mu0 = 0, V0 = 1)
  for (j in 1:4) {
    ess <- 150000 /
      (1 +
        2 * sum(stats::acf(g$chain[, j], plot = FALSE, lag.max = 100)$acf[-1]))
    z <- (mean(g$chain[, j]) - mean(g$fwd[, j])) /
      sqrt(var(g$fwd[, j]) / 150000 + var(g$chain[, j]) / ess)
    expect_lt(abs(z), 4.5)
  }
})

test_that("SV sampler recovers a known volatility path and coefficients; Sigma = A^-1 L A^-Tn", {
  set.seed(4)
  Tn <- 400
  A <- matrix(c(1, 0.5, 0, 1), 2, 2)
  Pi1 <- matrix(c(0.5, 0.1, 0, 0.4), 2, 2, byrow = TRUE)
  cc <- c(0.1, 0.2)
  h1 <- cumsum(rnorm(Tn, 0, 0.12))
  h1 <- h1 - mean(h1) - 0.5
  Y <- matrix(0, Tn, 2)
  for (t in 2:Tn) {
    Y[t, ] <- solve(
      A,
      cc + Pi1 %*% Y[t - 1, ] + exp(0.5 * c(h1[t], -0.7)) * rnorm(2)
    )
  }
  f <- bvarm(
    Y,
    p = 1,
    prior = prior_triangular("minnesota", lambda1 = 1, lambda2 = 1, delta = 0),
    sv = TRUE,
    ndraw = 1500,
    nburn = 1000
  )
  expect_gt(cor(f$h_mean[, 1], h1[-1]), 0.85)
  expect_lt(abs(mean(f$A[2, 1, ]) - 0.5), 0.1)
  expect_lt(
    max(abs(coef(f)[, 1] - solve(A, cbind(cc, Pi1))[1, ][c(2, 3, 1)])),
    0.25
  )
  d <- cpp_bvar_tri(
    cpp_var_ols(Y, 1, make_det(Tn, TRUE, FALSE, NULL))$Y,
    cpp_var_ols(Y, 1, make_det(Tn, TRUE, FALSE, NULL))$X,
    1L,
    1L,
    0L,
    TRUE,
    c(1, 1),
    c(0, 0),
    c(1, 1, 2, 100, 10, 0.1, 0.5),
    c(5, 0.04, 10),
    c(3, 2),
    3L,
    20L,
    1L,
    TRUE
  )
  for (k in 1:3) {
    Ai <- solve(d$A[,, k])
    expect_equal(
      d$Sigma[,, k],
      Ai %*% diag(exp(d$h_all[k, Tn - 1, ]), 2) %*% t(Ai),
      tolerance = 1e-12
    )
  }
})

test_that("SSVS and horseshoe separate true nonzero VAR coefficients", {
  set.seed(5)
  K <- 5
  Tn <- 250
  A1 <- diag(0.5, K)
  A1[cbind(c(2, 4), c(1, 3))] <- 0.3
  Y <- matrix(0, Tn + 50, K)
  for (t in 2:(Tn + 50)) {
    Y[t, ] <- A1 %*% Y[t - 1, ] + rnorm(K)
  }
  Y <- Y[-(1:50), ]
  nz <- c(t(A1)) != 0
  fs <- bvarm(
    Y,
    p = 2,
    prior = prior_triangular(
      "ssvs",
      lambda1 = 0.5,
      lambda2 = 1,
      alpha = 1,
      delta = 0,
      rho = 0.05
    ),
    ndraw = 1500,
    nburn = 800
  )
  inc <- fs$inclusion[1:K, ]
  expect_gt(mean(inc[matrix(nz, K)]), 0.85)
  expect_lt(mean(inc[!matrix(nz, K)]), 0.3)
  fh <- bvarm(
    Y,
    p = 2,
    prior = prior_triangular("horseshoe", lambda2 = 1, alpha = 1, delta = 0),
    ndraw = 1500,
    nburn = 800
  )
  b <- abs(coef(fh)[1:K, ])
  expect_lt(mean(b[!matrix(nz, K)]), 0.4 * mean(b[matrix(nz, K)]))
})

test_that("conjugate fixed-hyperparameter draws match analytic NIW moments", {
  set.seed(6)
  Y <- as.matrix(yraw[, 1:3])
  p <- 2
  M <- 3
  f <- bvarm(
    Y,
    p = p,
    prior = prior_minnesota(
      lambda = 0.3,
      soc = NULL,
      sur = NULL,
      hyper = character(0)
    ),
    ndraw = 30000
  )
  o <- cpp_var_ols(Y, p, make_det(nrow(Y), TRUE, FALSE, NULL))
  X <- o$X
  m <- ncol(X)
  psi <- ar_resid_var(Y, p)
  om <- c(0.3^2 / rep(1:p, each = M)^2 / rep(psi, p), 1e7)
  b0 <- matrix(0, m, M)
  b0[cbind(1:M, 1:M)] <- 1
  Gi <- solve(crossprod(X) + diag(1 / om))
  bh <- Gi %*% (crossprod(X, o$Y) + diag(1 / om) %*% b0)
  S <- diag(psi) +
    crossprod(o$Y - X %*% bh) +
    t(bh - b0) %*% diag(1 / om) %*% (bh - b0)
  nu <- nrow(o$Y) + M + 2
  expect_equal(
    apply(f$draws$Sigma, 1:2, mean),
    S / (nu - M - 1),
    tolerance = 0.02
  )
  expect_equal(
    unname(apply(f$draws$B, 1:2, mean)),
    unname(bh),
    tolerance = 0.05 / max(abs(bh))
  )
  expect_equal(
    diag(cov(t(matrix(f$draws$B, m * M)))),
    diag(kronecker(S / (nu - M - 1), Gi)),
    tolerance = 0.06
  )
})

test_that("Haar rotations are uniform and sign restrictions hold", {
  set.seed(7)
  N <- 40000
  drive <- function(K, S) {
    cpp_sign_irf(
      array(0, c(K, K, N)),
      array(diag(K), c(K, K, N)),
      1L,
      1L,
      S,
      1L,
      1L
    )
  }
  Q <- drive(2, array(0, c(2, 2, 1)))$P
  expect_lt(
    ks.test(Q[1, 1, ], function(x) {
      0.5 + asin(pmax(pmin(x, 1), -1)) / pi
    })$statistic,
    0.01
  )
  Q3 <- drive(3, array(0, c(3, 3, 1)))$P
  expect_lt(ks.test(Q3[1, 1, ], "punif", -1, 1)$statistic, 0.01)
  S2 <- array(0, c(2, 2, 1))
  S2[, 1, 1] <- c(1, 1)
  S2[, 2, 1] <- c(-1, 1)
  expect_equal(drive(2, S2)$accepted / N, 0.5, tolerance = 0.02)
  S3 <- array(0, c(3, 3, 1))
  S3[, 1, 1] <- 1
  expect_equal(drive(3, S3)$accepted / N, 0.25, tolerance = 0.03)
  R <- array(0, c(3, 3, 3))
  R[1, 1, ] <- 1
  R[2, 1, ] <- 1
  R[3, 2, ] <- -1
  fit <- bvarm(
    as.matrix(yraw[, 1:3]),
    p = 2,
    prior = prior_minnesota(hyper = character(0)),
    ndraw = 200
  )
  ir <- irf(fit, id = id_sign(R, nrot = 3, maxtry = 5000), horizon = 6)
  expect_true(
    all(ir$draws[1, 1, 1:3, ] > 0) &&
      all(ir$draws[2, 1, 1:3, ] > 0) &&
      all(ir$draws[3, 2, 1:3, ] < 0)
  )
})

test_that("Monte Carlo draws have the OLS asymptotic covariance", {
  set.seed(8)
  Y <- as.matrix(yraw[, 1:2])
  fit <- varm(Y, p = 2)
  d <- draws(fit, "mc", n = 20000)
  expect_equal(
    unname(apply(d$B, 1:2, mean)),
    unname(fit$B),
    tolerance = 0.02 / max(abs(fit$B))
  )
  expect_equal(
    unname(cov(t(matrix(d$B, length(fit$B))))[
      seq_len(nrow(fit$B)),
      seq_len(nrow(fit$B))
    ]),
    unname(fit$Sigma[1, 1] * fit$XtXinv),
    tolerance = 0.12
  )
})

test_that("mcmc_diagnostics: ESS, R-hat and Geweke behave", {
  set.seed(9)
  ess <- mean(replicate(
    15,
    mcmc_diagnostics(cbind(as.numeric(arima.sim(list(ar = 0.8), 3000))))$ess
  ))
  expect_equal(ess, 3000 * 0.2 / 1.8, tolerance = 0.12)
  ch <- lapply(1:4, function(j) cbind(a = rnorm(1000)))
  d <- mcmc_diagnostics(ch)
  expect_lt(abs(d$rhat - 1), 0.01)
  expect_gt(d$ess, 3000)
  expect_gt(
    mcmc_diagnostics(lapply(1:4, function(j) {
      cbind(a = rnorm(1000) + 3 * j)
    }))$rhat,
    1.5
  )
  expect_gt(
    abs(
      mcmc_diagnostics(cbind(cumsum(rnorm(2000, 0.02)) + rnorm(2000)))$geweke_z
    ),
    3
  )
  expect_s3_class(
    mcmc_diagnostics(bvarm(
      as.matrix(yraw[, 1:3]),
      p = 1,
      ndraw = 400,
      nburn = 200
    )),
    "data.frame"
  )
})

test_that("mcmc_diagnostics tidiers and plots", {
  set.seed(10)
  d <- mcmc_diagnostics(lapply(1:2, function(j) {
    cbind(a = rnorm(300), b = rnorm(300))
  }))
  expect_s3_class(d, "tbl_df")
  expect_named(d, c("term", "estimate", "std.error", "ess", "rhat", "geweke_z"))
  expect_equal(tidy(d), tibble::as_tibble(as.data.frame(d)))
  expect_equal(nrow(tidy(d)), 2L)
  gl <- glance(d)
  expect_equal(c(gl$n_param, gl$n_chain, gl$n_iter), c(2, 2, 300))
  for (type in c("summary", "trace", "density")) {
    p <- autoplot(d, type)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_s3_class(autoplot(d, "trace", params = b), "ggplot")
})
