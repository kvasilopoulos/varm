Yt <- as.data.frame(yraw)
colnames(Yt) <- c("a", "b", "c")
lagk <- function(v, k) c(rep(NA, k), v[seq_len(length(v) - k)])
leadk <- function(v, k) c(v[-seq_len(k)], rep(NA, k))

test_that("lp_joint: diagonal equals lp() se; cross-horizon covariance matches hand-coded HAC", {
  f <- lp_joint(
    Yt,
    "c",
    horizon = 3,
    lags = 2,
    contemp = c("a", "b"),
    responses = "a"
  )
  expect_equal(
    unname(sqrt(diag(f$vcov$a))),
    unname(f$se[, 1]),
    tolerance = 1e-10
  )
  expect_true(isSymmetric(unname(f$vcov$a)))
  # hand code: scores q_{t,h} e_{t+h}, aligned by t, h = 1, 2
  n <- nrow(Yt)
  d <- data.frame(
    Yt,
    a1 = lagk(Yt$a, 1),
    b1 = lagk(Yt$b, 1),
    c1 = lagk(Yt$c, 1),
    a2 = lagk(Yt$a, 2),
    b2 = lagk(Yt$b, 2),
    c2 = lagk(Yt$c, 2)
  )
  sc <- function(h) {
    dd <- d
    dd$y <- leadk(Yt$a, h)
    idx <- 3:(n - h)
    m <- lm(y ~ c + a1 + b1 + c1 + a2 + b2 + c2 + a + b, dd[idx, ])
    X <- model.matrix(m)
    q <- drop(X %*% solve(crossprod(X))[, 2]) * resid(m)
    out <- rep(0, n - 2)
    out[seq_along(q)] <- q
    out
  }
  g1 <- sc(1)
  g2 <- sc(2)
  L <- 3
  cr <- sum(g1 * g2) +
    sum(vapply(
      1:L,
      function(l) {
        (1 - l / (L + 1)) *
          (sum(g1[-(1:l)] * g2[1:(length(g2) - l)]) +
            sum(g2[-(1:l)] * g1[1:(length(g1) - l)]))
      },
      0
    ))
  expect_equal(unname(f$vcov$a[2, 3]), cr, tolerance = 1e-8)
})

test_that("lp_supt critical value lies between pointwise and Bonferroni", {
  f <- lp_joint(
    Yt,
    "c",
    horizon = 6,
    lags = 2,
    contemp = c("a", "b"),
    responses = "a"
  )
  b <- lp_supt(f, level = 0.9, nsim = 4000, seed = 1)
  cc <- attr(b, "crit")[[1]]
  expect_gt(cc, qnorm(0.95))
  expect_lt(cc, qnorm(1 - 0.05 / 6))
  expect_true(all(b$conf.low <= b$pw_low + 1e-12))
})

test_that("lp_state (contemporaneous) matches lm + Newey-West on the interacted regression", {
  skip_if_not_installed("sandwich")
  H <- 2
  st <- lp_state(
    Yt,
    "c",
    "a",
    H,
    lags = 2,
    contemp = c("a", "b"),
    gamma = 3,
    low_share = 0.3,
    responses = "b"
  )
  z <- as.numeric(scale(Yt$a))
  F1 <- lagk(1 / (1 + exp(3 * (z - quantile(z, 0.3)))), 1)
  n <- nrow(Yt)
  h <- 2
  idx <- 3:(n - h)
  W <- cbind(
    1,
    sapply(c("a", "b", "c"), function(v) lagk(Yt[[v]], 1)),
    sapply(c("a", "b", "c"), function(v) lagk(Yt[[v]], 2)),
    Yt$a,
    Yt$b
  )
  X <- cbind(Yt$c * (1 - F1), Yt$c * F1, W * (1 - F1), W * F1)[idx, ]
  m <- lm(leadk(Yt$b, h)[idx] ~ 0 + X)
  V <- sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE, adjust = FALSE)
  expect_equal(
    unname(st$estimate[h + 1, 1, "high"]),
    unname(coef(m)[1]),
    tolerance = 1e-8
  )
  expect_equal(
    unname(st$estimate[h + 1, 1, "low"]),
    unname(coef(m)[2]),
    tolerance = 1e-8
  )
  expect_equal(unname(st$se[h + 1, 1, "high"]), sqrt(V[1, 1]), tolerance = 1e-8)
  expect_equal(
    unname(st$se[h + 1, 1, "diff"]),
    sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2]),
    tolerance = 1e-8
  )
  expect_equal(
    st$estimate[,, "diff"],
    st$estimate[,, "low"] - st$estimate[,, "high"]
  )
})

test_that("lp_state (lagged timing) reproduces lpirfs::lp_nl point estimates", {
  skip_if_not_installed("lpirfs")
  H <- 4
  ref <- lpirfs::lp_nl(
    Yt,
    lags_endog_lin = 2,
    lags_endog_nl = 2,
    trend = 0,
    shock_type = 0,
    confint = 1,
    hor = H,
    switching = Yt$a,
    use_hp = 1,
    lambda = 1600,
    gamma = 3,
    use_nw = TRUE,
    num_cores = 1
  )
  me <- lp_state(
    Yt,
    2,
    "a",
    H,
    lags = 2,
    timing = "lagged",
    lags_lin = 2,
    shock_type = 0,
    gamma = 3,
    hp_lambda = 1600,
    state_const = FALSE
  )
  expect_equal(
    unname(me$estimate[,, "high"]),
    unname(t(ref$irf_s1_mean[,, 2])),
    tolerance = 1e-8
  )
  expect_equal(
    unname(me$estimate[,, "low"]),
    unname(t(ref$irf_s2_mean[,, 2])),
    tolerance = 1e-8
  )
})

test_that("lp_smooth: no penalty and a full basis reproduce LP; penalty shrinks variance", {
  H <- 5
  s0 <- lp_smooth(
    Yt,
    "c",
    "b",
    H,
    2,
    contemp = "a",
    nbasis = H + 1,
    lambda = 1e-10
  )
  expect_equal(unname(s0$estimate[, 1]), s0$lp, tolerance = 1e-6)
  f <- lp(Yt, "c", H, 2, contemp = "a", responses = "b")
  expect_equal(s0$lp, unname(f$estimate[, 1]), tolerance = 1e-8)
  expect_equal(s0$lp_se, unname(f$se[, 1]), tolerance = 1e-8)
  s1 <- lp_smooth(
    Yt,
    "c",
    "b",
    H,
    2,
    contemp = "a",
    lambda = 1e6 * s0$cv$lambda[1] + 1e3
  )
  expect_lt(s1$edf, s0$edf)
  s2 <- lp_smooth(Yt, "c", "b", H, 2, contemp = "a")
  expect_true(is.finite(s2$lambda) && all(s2$se > 0))
})

test_that("lp_panel matches lm with dummies and hand-coded cluster errors (unbalanced)", {
  set.seed(2)
  N <- 10
  Tn <- 20
  d <- expand.grid(id = 1:N, t = 1:Tn)
  d$s <- rnorm(nrow(d))
  a <- rnorm(N)
  b <- rnorm(Tn)
  d$y <- 0.4 * d$s + a[d$id] + b[d$t] + rnorm(nrow(d))
  d <- d[order(d$id, d$t), ]
  d <- d[-sample(nrow(d), 15), ]
  res <- lp_panel(d, "y", "s", "id", "t", horizon = 2, lags = 1)
  key <- d$id * 1000 + d$t
  lead <- function(v, k) v[match(key + k, key)]
  for (h in 0:2) {
    dd <- d
    dd$yh <- lead(d$y, h)
    dd$yl <- lead(d$y, -1)
    dd$sl <- lead(d$s, -1)
    dd <- dd[complete.cases(dd), ]
    m <- lm(yh ~ s + yl + sl + factor(id) + factor(t), dd)
    X <- model.matrix(m)
    n <- nrow(X)
    k <- ncol(X)
    cs <- rowsum(X * resid(m), dd$id)
    G <- nrow(cs)
    Ri <- solve(crossprod(X))
    V <- G / (G - 1) * (n - 1) / (n - k) * Ri %*% crossprod(cs) %*% Ri
    expect_equal(res$estimate[h + 1, 1], unname(coef(m)["s"]), tolerance = 1e-8)
    expect_equal(res$se[h + 1, 1], sqrt(V[2, 2]), tolerance = 1e-8)
  }
  expect_equal(
    unname(sqrt(diag(res$vcov$y))),
    unname(res$se[, 1]),
    tolerance = 1e-10
  )
})

test_that("lp_forecast equals the VAR forecast at h = 1 and compare_lp_var is consistent", {
  fc <- lp_forecast(Yt, 3, 2)
  v <- forecast(varm(Yt, 2), h = 1, ci = "none")
  expect_equal(unname(fc$mean[1, ]), unname(v$mean[1, ]), tolerance = 1e-8)
  cmp <- compare_lp_var(Yt, "b", horizon = 5, lags = 2)
  expect_equal(cmp$diff, cmp$lp - cmp$var)
  expect_equal(cmp$lp[cmp$horizon == 0 & cmp$response == "b"], 1)
  expect_equal(cmp$var[cmp$horizon == 0 & cmp$response == "a"], 0)
})

test_that("tidiers of lp_joint, lp_supt, lp_state, lp_smooth, lp_panel, lp_forecast", {
  fj <- lp_joint(Yt, "c", horizon = 3, lags = 2, contemp = c(a, b))
  expect_equal(tidy(fj)$estimate, tidy.lpm(fj)$estimate)
  expect_equal(nrow(tidy(fj)), 4 * 3 * 2)
  sp <- lp_supt(fj, level = 0.9, nsim = 500, seed = 1)
  expect_s3_class(sp, "tbl_df")
  expect_named(
    sp,
    c(
      "response",
      "horizon",
      "estimate",
      "std.error",
      "pw_low",
      "pw_high",
      "conf.low",
      "conf.high"
    )
  )
  ts <- tidy(sp)
  expect_equal(nrow(ts), 2 * nrow(sp))
  expect_setequal(as.character(ts$band), c("pointwise", "sup-t"))

  st <- lp_state(
    Yt,
    c,
    state = a,
    horizon = 3,
    lags = 2,
    contemp = c(a, b),
    gamma = 3,
    low_share = 0.3
  )
  t1 <- tidy(st)
  expect_named(
    t1,
    c(
      "horizon",
      "response",
      "state",
      "estimate",
      "std.error",
      "level",
      "conf.low",
      "conf.high"
    )
  )
  expect_equal(levels(t1$state), c("high", "low", "diff"))
  expect_equal(nrow(t1), 4 * 3 * 3 * 2)
  expect_equal(nrow(glance(st)), 1L)
  expect_true("transition" %in% names(glance(st)))

  sm <- lp_smooth(Yt, c, b, horizon = 4, lags = 2, contemp = a)
  t2 <- tidy(sm, raw = TRUE)
  expect_setequal(t2$type, c("smooth", "raw"))
  expect_equal(sum(t2$type == "raw"), 5L)
  expect_true(all(c("lambda", "edf") %in% names(glance(sm))))

  set.seed(3)
  d <- expand.grid(id = 1:8, t = 1:25)
  d$s <- rnorm(nrow(d))
  d$y <- 0.5 * d$s + rnorm(nrow(d))
  pn <- lp_panel(d, y, s, id, t, horizon = 2, lags = 1)
  expect_equal(nrow(tidy(pn)), 3 * 2)
  expect_true(all(c("fe", "cluster", "nclusters") %in% names(glance(pn))))
  expect_equal(
    pn$estimate,
    lp_panel(d, "y", "s", "id", "t", horizon = 2, lags = 1)$estimate
  )

  fc <- tidy(lp_forecast(Yt, 3, 2))
  expect_named(
    fc,
    c(
      "horizon",
      "variable",
      "estimate",
      "std.error",
      "level",
      "conf.low",
      "conf.high"
    )
  )
  expect_equal(nrow(fc), 3 * 3 * 2)
})

test_that("autoplot methods for the lp family return ggplots that build", {
  fj <- lp_joint(Yt, "c", horizon = 3, lags = 2, contemp = c(a, b))
  st <- lp_state(
    Yt,
    c,
    state = a,
    horizon = 3,
    lags = 2,
    contemp = c(a, b),
    gamma = 3,
    low_share = 0.3
  )
  sm <- lp_smooth(Yt, c, b, horizon = 4, lags = 2, contemp = a)
  plots <- list(
    autoplot(fj),
    autoplot(st),
    autoplot(st, "diff"),
    autoplot(sm),
    autoplot(lp_supt(fj, nsim = 500, seed = 1)),
    autoplot(lp_forecast(Yt, 3, 2))
  )
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_s3_class(plot(st), "ggplot")
})

test_that("tidyselect works for state, smooth and panel arguments", {
  s1 <- lp_state(Yt, "c", "a", 2, lags = 2, responses = "b")
  s2 <- lp_state(Yt, c, a, 2, lags = 2, responses = starts_with("b"))
  expect_equal(s1$estimate, s2$estimate)
  expect_error(lp_state(Yt, c, zz, 2, lags = 2), "zz")
  expect_error(lp_smooth(Yt, c, c(a, b), 2, lags = 1), "exactly one")
  expect_error(compare_lp_var(Yt, zz), "zz")
})
