test_that("hierarchical NIW log-posterior matches BVAR at its optimum", {
  skip_if_not_installed("BVAR")
  p <- 4
  set.seed(1)
  bv <- suppressWarnings(BVAR::bvar(
    Y3,
    lags = p,
    n_draw = 400,
    n_burn = 100,
    n_thin = 1,
    priors = BVAR::bv_priors(
      hyper = "auto",
      mn = BVAR::bv_mn(),
      soc = BVAR::bv_soc(),
      sur = BVAR::bv_sur()
    ),
    verbose = FALSE
  ))
  psi_b <- bv$priors$psi$mode
  ols <- cpp_var_ols(Y3, p, matrix(1, nrow(Y3), 1))
  b <- matrix(0, 13, 3)
  b[cbind(1:3, 1:3)] <- 1
  gp <- rbind(
    gamma_coef(0.2, 0.4),
    gamma_coef(2, 0.25),
    gamma_coef(1, 1),
    gamma_coef(1, 1)
  )
  h <- c(
    bv$optim$par[["lambda"]],
    2,
    bv$optim$par[["soc"]],
    bv$optim$par[["sur"]],
    psi_b
  )
  v <- cpp_niw_logpost(
    h,
    ols$Y,
    ols$X,
    p,
    1L,
    TRUE,
    b,
    1e7,
    c(1L, 0L, 1L, 1L, 0L, 0L, 0L),
    gp,
    c(1e-4, 1, 1e-4, 1e-4, 0, 0, 0),
    c(5, 3, 50, 50, Inf, Inf, Inf),
    c(0.004, 0.004),
    FALSE
  )
  expect_equal(v, bv$optim$value, tolerance = 1e-6)
})

test_that("bvarm: fixed, hierarchical and triangular samplers produce sane draws", {
  set.seed(1)
  ff <- bvarm(Y3, 2, prior_minnesota(hyper = character(0)), ndraw = 200)
  expect_equal(dim(ff$draws$B), c(7, 3, 200))
  fh <- bvarm(Y3, 2, ndraw = 300, nburn = 150)
  expect_equal(colnames(fh$hyper), c("lambda", "soc", "sur"))
  expect_true(fh$accept > 0.05 && fh$accept < 0.95)
  expect_true(max(stability(fh)[1, ]) < 1.1)
  expect_output(print(fh), "bvarm")
  for (cfg in list(
    list("minnesota", FALSE),
    list("ssvs", FALSE),
    list("horseshoe", TRUE)
  )) {
    ft <- bvarm(
      Y3,
      2,
      prior_triangular(cfg[[1]]),
      ndraw = 200,
      nburn = 100,
      sv = cfg[[2]]
    )
    expect_equal(dim(ft$draws$B), c(7, 3, 200))
    expect_true(all(is.finite(coef(ft))))
    expect_s3_class(irf(ft, horizon = 4), "varm_irf")
  }
  expect_equal(dim(volatility(ft)), c(213, 3))
  expect_error(volatility(ff), "stochastic volatility")
  expect_error(bvarm(Y3, 2, sv = TRUE), "prior_triangular")
})

test_that("bvarm methods: irf, fevd, hd, forecast, K = 1 and exog", {
  set.seed(1)
  fh <- bvarm(Y3, 2, ndraw = 150, nburn = 100)
  expect_equal(dim(irf(fh, horizon = 3)$bands), c(3, 3, 4, 4))
  expect_s3_class(fevd(fh, horizon = 3), "varm_fevd")
  expect_s3_class(hd(fh, ndraw = 5), "varm_hd")
  fc <- forecast(fh, h = 3)
  expect_equal(dim(fc$bands), c(3, 3, 4))
  b1 <- bvarm(Y3[, 1, drop = FALSE], 2, ndraw = 100, nburn = 50)
  expect_s3_class(hd(b1, ndraw = 5), "varm_hd")
  expect_s3_class(irf(b1, horizon = 2), "varm_irf")
  be <- bvarm(Y3, 2, ndraw = 100, nburn = 50, exog = cbind(z = rnorm(nrow(Y3))))
  expect_error(forecast(be, h = 2), "exog_future")
  expect_s3_class(
    forecast(be, h = 2, exog_future = matrix(0, 2, 1)),
    "varm_forecast"
  )
})

test_that("hierarchical sampler works with soc/sur dropped (bounds of disabled hyperparameters)", {
  set.seed(1)
  f <- bvarm(
    yraw,
    2,
    prior_minnesota(soc = NULL, sur = NULL, hyper = "lambda"),
    ndraw = 300,
    nburn = 150
  )
  expect_s3_class(f, "bvarm")
  expect_equal(colnames(f$hyper), "lambda")
  expect_true(all(is.finite(f$hyper)))
  expect_gt(f$accept, 0.05)
})

test_that("tidy, glance and autoplot for bvarm", {
  set.seed(2)
  fh <- bvarm(Y3, 2, ndraw = 200, nburn = 100)
  td <- tidy(fh)
  expect_s3_class(td, "tbl_df")
  expect_named(
    td,
    c("response", "term", "estimate", "std.error", "conf.low", "conf.high")
  )
  expect_equal(nrow(td), 7 * 3)
  expect_equal(td$estimate, as.vector(coef(fh)))
  expect_true(all(td$conf.low <= td$estimate & td$estimate <= td$conf.high))
  expect_true(all(
    tidy(fh, conf.level = 0.5)$conf.high -
      tidy(fh, conf.level = 0.5)$conf.low <=
      td$conf.high - td$conf.low + 1e-12
  ))
  gl <- glance(fh)
  expect_equal(nrow(gl), 1L)
  expect_true(all(
    c(
      "K",
      "p",
      "nobs",
      "ndraw",
      "method",
      "logml",
      "accept",
      "lambda",
      "soc",
      "sur"
    ) %in%
      names(gl)
  ))
  expect_equal(gl$ndraw, 200L)
  for (type in c("coef", "hyper", "trace")) {
    p <- autoplot(fh, type)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  ff <- bvarm(Y3, 2, prior_minnesota(hyper = character(0)), ndraw = 100)
  expect_true(is.na(glance(ff)$accept))
  expect_error(autoplot(ff, "hyper"), "hyperparameter")
  expect_s3_class(autoplot(ff, "trace"), "ggplot")
  expect_s3_class(plot(fh), "ggplot")
})

test_that("bvarm accepts data frames with extra columns and tidyselect", {
  set.seed(3)
  df <- tibble::as_tibble(yraw)
  df$date <- as.Date("2000-01-01") + seq_len(nrow(df))
  f <- bvarm(
    df,
    1,
    prior_minnesota(hyper = character(0)),
    ndraw = 50,
    vars = c(infl, rate)
  )
  expect_equal(f$varnames, c("infl", "rate"))
})
