test_that("tvpvar runs and has consistent dimensions", {
  set.seed(1)
  f <- tvpvar(yraw, p = 2, ndraw = 100, nburn = 100)
  expect_s3_class(f, "tvpvar")
  expect_equal(dim(f$beta)[1:3], c(7, 3, 173))
  expect_equal(dim(f$beta)[4], 100)
  expect_equal(dim(f$h), c(3, 173, 100))
  expect_true(all(is.finite(f$beta)) && all(is.finite(f$h)))
  expect_output(print(f), "TVP-VAR")
  expect_equal(dim(coef(f, date = "last")), c(7, 3))
  expect_s3_class(plot(f, "coef"), "ggplot")
  expect_s3_class(plot(f, "vol"), "ggplot")
})

test_that("constant limit recovers OLS", {
  set.seed(2)
  pr <- tvp_prior(tau = 40, kQ = 1e-6, kS = 1e-6, kW = 1e-6, vscale = 1e4)
  f <- tvpvar(
    yraw,
    p = 2,
    ndraw = 1500,
    nburn = 500,
    sv = FALSE,
    tv_A = FALSE,
    prior = pr
  )
  ols <- varm(yraw[41:215, ], p = 2)
  expect_lt(max(abs(coef(f, date = "last") - ols$B)), 0.03)
  expect_lt(max(apply(f$beta_mean, 1:2, function(v) diff(range(v)))), 1e-4)
  fc <- tvpvar(
    yraw,
    p = 2,
    ndraw = 1500,
    nburn = 500,
    sv = FALSE,
    tv_beta = FALSE,
    tv_A = FALSE,
    prior = pr
  )
  expect_lt(max(abs(coef(fc, date = "last") - ols$B)), 0.03)
})

test_that("irf, fevd, forecast and paths", {
  set.seed(3)
  f <- tvpvar(yraw, p = 2, ndraw = 100, nburn = 100)
  i <- irf(f, date = c(60, 200), horizon = 6)
  expect_s3_class(i, "tvpvar_irf")
  expect_length(i, 2)
  i1 <- irf(f, date = "last", horizon = 6)
  expect_s3_class(i1, "varm_irf")
  expect_equal(dim(i1$estimate), c(3, 3, 7))
  expect_true(all(i1$estimate[,, 1][upper.tri(diag(3))] == 0)) # recursive impact
  in1 <- irf(
    f,
    date = "last",
    horizon = 4,
    normalise = list(shock = 3, variable = 3, size = 1)
  )
  expect_equal(unname(in1$estimate[3, 3, 1]), 1)
  fv <- fevd(f, date = "last", horizon = 6)
  expect_equal(unname(rowSums(fv$estimate[,, 7])), rep(1, 3), tolerance = 0.05)
  fc <- forecast(f, h = 5)
  expect_equal(dim(fc$draws), c(100, 5, 3))
  expect_true(all(is.finite(fc$draws)))
  ip <- irf_path(f, shock = 3, response = 3, horizon = 0)
  expect_equal(nrow(ip), 173)
  expect_named(ip, c("date", "estimate", "conf.low", "conf.high"))
  expect_true(all(ip$conf.high >= ip$conf.low))
  ip2 <- irf_path(f, shock = rate, response = rate, horizon = 0)
  expect_equal(ip2$estimate, ip$estimate)
  expect_s3_class(plot(ip), "ggplot")
  vol <- tvp_volatility(f)
  expect_equal(nrow(vol), 3 * 173)
  expect_named(vol, c("date", "variable", "estimate", "conf.low", "conf.high"))
})

test_that("variants run (Cogley-Sargent, diagonal Q, homoskedastic)", {
  set.seed(4)
  for (a in list(
    list(tv_A = FALSE),
    list(sv = FALSE, tv_A = FALSE),
    list(q_diag = TRUE)
  )) {
    f <- do.call(tvpvar, c(list(yraw, p = 1, ndraw = 50, nburn = 50), a))
    expect_true(all(is.finite(f$beta)))
  }
  f <- tvpvar(yraw, p = 1, ndraw = 30, nburn = 30, reject_explosive = TRUE)
  expect_true(all(is.finite(f$beta)))
})

test_that("sampler is calibrated when the DGP is the prior (K=2)", {
  set.seed(5)
  K <- 2
  m <- 3
  nb <- 6
  Tn <- 100
  N <- 200
  b0 <- c(0.5, 0.1, 0, 0.1, 0.5, 0)
  P0 <- diag(0.01, nb)
  Qv <- 1e-4
  Sv <- 1e-3
  Wv <- 0.01
  hit <- c(0, 0, 0)
  for (r in 1:6) {
    b <- b0 + 0.1 * rnorm(nb)
    a <- 0.3 + 0.2 * rnorm(1)
    h <- log(0.25) + sqrt(0.5) * rnorm(2)
    Y <- matrix(0, Tn, K)
    X <- matrix(0, Tn, m)
    y <- c(0, 0)
    TB <- matrix(0, nb, Tn)
    TA <- numeric(Tn)
    TH <- matrix(0, 2, Tn)
    for (t in 1:Tn) {
      b <- b + sqrt(Qv) * rnorm(nb)
      a <- a + sqrt(Sv) * rnorm(1)
      h <- h + sqrt(Wv) * rnorm(2)
      x <- c(y, 1)
      X[t, ] <- x
      y <- as.numeric(
        x %*%
          matrix(b, m, K) +
          solve(
            matrix(c(1, 0, a, 1), 2, 2, byrow = TRUE),
            exp(h / 2) * rnorm(2)
          )
      )
      Y[t, ] <- y
      TB[, t] <- b
      TA[t] <- a
      TH[, t] <- h
    }
    o <- varm:::cpp_tvp(
      Y,
      X,
      1L,
      b0,
      P0,
      0.3,
      list(NULL, matrix(0.04)),
      diag(Qv * 1e6, nb),
      1e6,
      diag(Qv, nb),
      list(NULL, matrix(Sv * 1e6)),
      c(0, 1e6),
      list(NULL, matrix(Sv)),
      rep(log(0.25), 2),
      0.5,
      c(1e5, 1e5 * Wv),
      c(3, 2),
      c(1, 1, 1, 0, 0, 1),
      N,
      200L,
      1L,
      N
    )
    cv <- function(a, tr) {
      q <- apply(a, 1:2, quantile, probs = c(0.05, 0.95))
      mean(tr >= q[1, , ] & tr <= q[2, , ])
    }
    hit <- hit +
      c(
        cv(array(o$beta, c(nb, Tn, N)), TB),
        cv(array(o$alpha, c(1, Tn, N)), matrix(TA, 1)),
        cv(o$h, TH)
      )
  }
  expect_true(all(abs(hit / 6 - 0.9) < 0.12))
})

test_that("tidy, glance and autoplot for tvpvar", {
  set.seed(6)
  f <- tvpvar(yraw, p = 1, ndraw = 60, nburn = 60)
  tc <- tidy(f)
  expect_s3_class(tc, "tbl_df")
  expect_named(
    tc,
    c("date", "response", "term", "estimate", "conf.low", "conf.high")
  )
  expect_equal(nrow(tc), f$T * 3 * 4)
  expect_true(all(tc$conf.low <= tc$conf.high))
  expect_equal(nrow(tidy(f, equation = infl, const = FALSE)), f$T * 3)
  expect_named(
    tidy(f, "vol"),
    c("date", "variable", "estimate", "conf.low", "conf.high")
  )
  expect_equal(tibble::as_tibble(f), tc)
  gl <- glance(f)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$nobs, f$T)
  for (p in list(
    autoplot(f),
    autoplot(f, "coef", equation = c(infl, rate)),
    autoplot(f, "vol"),
    autoplot(f, "shock_vol"),
    autoplot(f, "irf", shock = rate, response = infl, horizon = c(0, 4)),
    autoplot(irf_path(f, rate, infl, horizon = 2))
  )) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_error(autoplot(f, "irf"), "shock")
})

test_that("tvpvar accepts tidyselect variables and a Date column", {
  set.seed(7)
  df <- tibble::as_tibble(yraw)
  df$date <- seq(as.Date("1960-01-01"), by = "quarter", length.out = nrow(df))
  f <- tvpvar(df, p = 1, ndraw = 40, nburn = 40, vars = c(infl, rate))
  expect_equal(f$varnames, c("infl", "rate"))
  expect_s3_class(f$dates, "Date")
  expect_equal(nrow(tidy(f, "vol")), 2 * f$T)
  expect_s3_class(autoplot(f, "vol"), "ggplot")
})

test_that("tvpvar warns when the Q prior has undefined mean and irf() rejects `id`", {
  expect_warning(
    f4 <- tvpvar(yraw, p = 4, ndraw = 30, nburn = 10),
    "prior degrees of freedom for Q"
  )
  expect_no_warning(
    f2 <- tvpvar(yraw, p = 2, ndraw = 30, nburn = 10)
  )
  expect_error(irf(f2, id = id_chol(), horizon = 4), "recursive")
  expect_warning(
    irf(f2, date = 1800, horizon = 4),
    "outside the estimation sample"
  )
})
