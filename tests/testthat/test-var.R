skip_if_not_installed("vars")
f <- varm(Y3, p = 2)
v <- vars::VAR(Y3, p = 2, type = "const")

test_that("OLS, logLik, IC and stability agree with vars", {
  Bv <- sapply(v$varresult, coef)[c(1:6, 7), ]
  expect_equal(unname(f$B), unname(Bv), tolerance = 1e-9)
  expect_equal(unname(f$Sigma), unname(summary(v)$covres), tolerance = 1e-9)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(v)), tolerance = 1e-9)
  vs <- vars::VARselect(Y3, lag.max = 6, type = "const")
  mine <- var_select(Y3, 6)
  expect_equal(unname(mine$selection), unname(vs$selection))
  expect_equal(
    unname(mine$criteria[, 1:3]),
    unname(t(vs$criteria)[, 1:3]),
    tolerance = 1e-9
  )
  expect_equal(
    sort(stability(f), decreasing = TRUE),
    sort(Mod(vars::roots(v)), decreasing = TRUE),
    tolerance = 1e-9
  )
})

test_that("Cholesky IRF and FEVD agree with vars", {
  ir <- vars::irf(v, n.ahead = 12, ortho = TRUE, boot = FALSE)$irf
  m <- irf(f, horizon = 12, ci = "none")$estimate
  fe <- vars::fevd(v, n.ahead = 13)
  mf <- fevd(f, horizon = 12, ci = "none")$estimate
  for (s in colnames(Y3)) {
    for (r in colnames(Y3)) {
      expect_equal(unname(m[r, s, ]), unname(ir[[s]][, r]), tolerance = 1e-9)
      expect_equal(unname(mf[r, s, ]), unname(fe[[r]][, s]), tolerance = 1e-9)
    }
  }
})

test_that("Blanchard-Quah impact matrix agrees with vars (up to sign)", {
  bq <- vars::BQ(vars::VAR(Y3[, 1:2], p = 2))
  P <- irf(
    varm(Y3[, 1:2], p = 2),
    id = id_longrun(),
    horizon = 0,
    ci = "none"
  )$estimate[,, 1]
  expect_equal(unname(abs(P)), unname(abs(bq$B)), tolerance = 1e-8)
})

test_that("point forecast agrees with vars; exog forecast follows the recursion", {
  fc <- predict(v, n.ahead = 8, ci = 0.95)$fcst
  mf <- forecast(f, h = 8, ci = "none")$mean
  for (k in colnames(Y3)) {
    expect_equal(unname(mf[, k]), unname(fc[[k]][, "fcst"]), tolerance = 1e-9)
  }
  set.seed(1)
  X <- cbind(z = rnorm(nrow(Y3)))
  fe <- varm(Y3, p = 2, exog = X, trend = TRUE)
  xf <- cbind(z = c(1, -1, 0.5))
  m <- forecast(fe, h = 3, exog_future = xf, ci = "none")$mean
  Tn <- nrow(Y3)
  x1 <- c(Y3[Tn, ], Y3[Tn - 1, ], 1, Tn + 1, xf[1, 1])
  expect_equal(unname(m[1, ]), unname(drop(x1 %*% fe$B)), tolerance = 1e-10)
  expect_error(forecast(fe, h = 3), "exog_future")
  expect_error(
    forecast(fe, h = 3, exog_future = xf[1:2, , drop = FALSE]),
    "h x q"
  )
})

test_that("horizon 0, normalise, cumulative and sign restrictions with bands", {
  i0 <- irf(f, horizon = 0, ci = "none")
  expect_equal(dim(i0$estimate), c(3, 3, 1))
  i0b <- irf(f, horizon = 0, ci = "mc", nboot = 30)
  expect_equal(dim(i0b$bands), c(3, 3, 1, 4))
  expect_equal(dim(fevd(f, horizon = 0, ci = "none")$estimate), c(3, 3, 1))
  n <- irf(
    f,
    horizon = 4,
    ci = "mc",
    nboot = 30,
    normalise = list(shock = "rate", variable = "rate", size = 2)
  )
  expect_equal(unname(n$estimate["rate", "rate", 1]), 2)
  cu <- irf(f, horizon = 4, ci = "none", cumulative = "infl")
  ref <- irf(f, horizon = 4, ci = "none")
  expect_equal(
    unname(cu$estimate["infl", , ]),
    unname(t(apply(ref$estimate["infl", , ], 1, cumsum)))
  )
  S <- matrix(c(1, 1, 1, 1, -1, 0, 1, 0, 0), 3, 3, byrow = TRUE)
  set.seed(1)
  sg <- irf(f, id = id_sign(S), horizon = 4, ci = "mc", nboot = 30)
  expect_false(is.null(sg$bands))
  expect_equal(dim(sg$set_estimate), c(3, 3, 5))
  set.seed(1)
  expect_s3_class(
    fevd(f, id = id_sign(S), horizon = 4, ci = "mc", nboot = 20),
    "varm_fevd"
  )
  expect_error(irf(f, id = id_sign(diag(2)), ci = "none"), "dimension")
})

test_that("bootstrap variants, including bias correction, run and bracket the estimate", {
  set.seed(1)
  for (ci in c("bootstrap", "wild", "wild_mammen")) {
    r <- irf(f, horizon = 4, ci = ci, nboot = 60)
    expect_equal(dim(r$bands), c(3, 3, 5, 4))
    expect_true(all(is.finite(r$bands)))
  }
  r <- irf(f, horizon = 4, ci = "bootstrap", nboot = 60, bias_correct = TRUE)
  expect_true(all(is.finite(r$bands)))
})

test_that("edge cases: K = 1, p = 1, deterministic options, ts/data.frame input", {
  f1 <- varm(Y3[, 1], p = 2)
  expect_equal(f1$K, 1L)
  expect_equal(
    dim(irf(f1, horizon = 3, ci = "mc", nboot = 20)$estimate),
    c(1, 1, 4)
  )
  expect_equal(dim(fevd(f1, horizon = 3, ci = "none")$estimate), c(1, 1, 4))
  expect_equal(dim(hd(f1)$contributions)[2:3], c(1, 1))
  expect_equal(nrow(forecast(f1, h = 3, ci = "none")$mean), 3)
  f2 <- varm(Y3, p = 1)
  expect_equal(
    dim(irf(f2, horizon = 3, ci = "bootstrap", nboot = 20)$estimate),
    c(3, 3, 4)
  )
  expect_equal(nrow(forecast(f2, h = 1, ci = "mc", nboot = 20)$mean), 1)
  fn <- varm(Y3, p = 2, const = FALSE)
  expect_equal(rownames(fn$B), rownames(f$B)[1:6])
  expect_equal(nrow(forecast(fn, h = 2, ci = "none")$mean), 2)
  ft <- varm(Y3, p = 2, trend = TRUE)
  expect_true("trend" %in% rownames(ft$B))
  yts <- ts(Y3, start = c(1990, 1), frequency = 4)
  expect_equal(unname(varm(yts, 2)$B), unname(f$B))
  expect_equal(unname(varm(yraw, 2)$B), unname(f$B))
})

test_that("bad input gives informative errors", {
  expect_error(varm(Y3, p = 0), "`p` must be >= 1")
  expect_error(varm(Y3, p = 300), "not enough observations")
  expect_error(varm("a"), "numeric")
  Yna <- Y3
  Yna[3, 2] <- NA
  expect_error(varm(Yna, 2), "missing")
  expect_error(varm(Y3, 2, exog = matrix(1, 3, 1)), "same number of rows")
  expect_error(irf(f, id = id_chol("zzz"), ci = "none"), "unknown")
  expect_error(irf(f, id = id_chol(zzz), ci = "none"), "unknown")
  expect_error(
    irf(f, id = id_chol(c("rate", "infl")), ci = "none"),
    "permutation"
  )
  expect_error(
    irf(f, normalise = list(shock = "zzz", variable = "rate"), ci = "none"),
    "unknown shock"
  )
  expect_error(irf(f, cumulative = "zzz", ci = "none"), "unknown")
  expect_error(id_sign(matrix(2, 3, 3)), "-1, 0 or \\+1")
  expect_error(irf(f, ci = "foo"))
})

test_that("tidyselect works in vars, order, cumulative and normalise_to", {
  fv <- varm(yraw, p = 2, vars = c(infl, rate))
  expect_equal(fv$varnames, c("infl", "rate"))
  expect_equal(varm(yraw, p = 2, vars = c("infl", "rate"))$B, fv$B)
  expect_equal(varm(yraw, p = 2, vars = 1:2)$varnames, c("infl", "unemp"))
  expect_equal(
    var_select(yraw, 3, vars = starts_with("i"))$criteria[, 1],
    var_select(yraw[1], 3)$criteria[, 1]
  )
  a <- irf(
    f,
    id = id_chol(order = c(rate, infl, unemp)),
    horizon = 3,
    ci = "none"
  )
  b <- irf(
    f,
    id = id_chol(c("rate", "infl", "unemp")),
    horizon = 3,
    ci = "none"
  )
  expect_equal(a$estimate, b$estimate)
  ref <- irf(f, horizon = 4, ci = "none")
  cu <- irf(f, horizon = 4, ci = "none", cumulative = c(infl, rate))
  expect_equal(
    unname(cu$estimate["rate", , ]),
    unname(t(apply(ref$estimate["rate", , ], 1, cumsum)))
  )
  n1 <- irf(
    f,
    horizon = 2,
    ci = "none",
    normalise = normalise_to(rate, rate, 2)
  )
  n2 <- irf(
    f,
    horizon = 2,
    ci = "none",
    normalise = list(shock = "rate", variable = "rate", size = 2)
  )
  expect_equal(n1$estimate, n2$estimate)
  expect_error(
    irf(f, ci = "none", normalise = normalise_to(rate, zzz)),
    "unknown variable"
  )
  expect_error(
    irf(f, ci = "none", normalise = normalise_to(c(rate, infl), rate)),
    "exactly one"
  )
})

test_that("data frames: non-numeric columns dropped, Date column kept", {
  d <- cbind(
    date = seq(as.Date("1990-01-01"), by = "quarter", length.out = nrow(yraw)),
    yraw
  )
  fd <- varm(d, p = 2)
  expect_equal(unname(fd$B), unname(f$B))
  expect_s3_class(fd$dates, "Date")
  expect_equal(fd$varnames, names(yraw))
  expect_error(varm(data.frame(a = letters[1:20]), 1), "no numeric columns")
})
