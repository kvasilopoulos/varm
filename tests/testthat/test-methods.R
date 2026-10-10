f <- varm(Y3, p = 2)

test_that("print/summary/accessors work", {
  expect_output(print(f), "VAR\\(2\\)")
  expect_output(print(summary(f)), "Equation: infl")
  expect_equal(dim(coef(f)), c(7, 3))
  expect_equal(dim(resid(f)), c(213, 3))
  expect_equal(dim(fitted(f)), c(213, 3))
  expect_equal(dim(vcov(f)), c(3, 3))
  expect_equal(nobs(f), 213L)
  expect_s3_class(logLik(f), "logLik")
  expect_equal(dim(lag_matrices(f)), c(3, 3, 2))
  expect_output(print(draws(f, "mc", 20)), "varm_draws")
  expect_equal(dim(stability(draws(f, "mc", 20))), c(20, 6))
})

test_that("IRF object methods: print, as.data.frame, plot, autoplot", {
  r <- irf(f, horizon = 5, ci = "mc", nboot = 30)
  expect_output(print(r), "IRF")
  df <- as.data.frame(r)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 9 * 6 * 2)
  expect_true(all(c("response", "shock", "horizon", "estimate") %in% names(df)))
  expect_setequal(unique(df$level), c(0.68, 0.9))
  expect_equal(nrow(as.data.frame(r, bands = FALSE)), 9 * 6)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plot(r), "ggplot")
  expect_s3_class(
    plot(r, shocks = "rate", responses = c("infl", "unemp")),
    "ggplot"
  )
  expect_s3_class(ggplot2::autoplot(r), "ggplot")
  expect_s3_class(plot(fevd(f, horizon = 5, ci = "none")), "ggplot")
  expect_output(print(id_chol()), "chol", ignore.case = TRUE)
  expect_output(print(id_chol(order = c(rate, infl, unemp))), "rate, infl")
})

test_that("hd, forecast and lpm methods", {
  h <- hd(f)
  expect_equal(dim(h$contributions), c(213, 3, 3))
  expect_equal(
    unname(apply(h$contributions, c(1, 2), sum) + h$baseline),
    unname(h$actual),
    tolerance = 1e-8
  )
  fc <- forecast(f, h = 4)
  expect_output(print(fc), "forecast")
  expect_equal(dim(fc$bands), c(4, 3, 4))
  l <- lp(Y3, "rate", horizon = 8, lags = 2)
  expect_output(print(l), "lpm")
  expect_s3_class(as.data.frame(l), "data.frame")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(l))
})

test_that("yraw has the documented shape", {
  expect_named(yraw, c("infl", "unemp", "rate"))
  expect_s3_class(yraw, "data.frame")
})
