fit <- varm(yraw, p = 2)

test_that("tidy/glance/augment for varm follow broom conventions", {
  td <- tidy(fit)
  expect_s3_class(td, "tbl_df")
  expect_named(
    td,
    c("response", "term", "estimate", "std.error", "statistic", "p.value")
  )
  expect_length(fit$B, nrow(td))
  s <- summary(fit)
  expect_equal(td$estimate, as.vector(fit$B))
  expect_equal(td$p.value, as.vector(s$pvalue))
  ci <- tidy(fit, conf.int = TRUE, conf.level = 0.9)
  expect_true(all(ci$conf.low < ci$estimate & ci$estimate < ci$conf.high))
  g <- glance(fit)
  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1)
  expect_named(g, c("K", "p", "nobs", "logLik", "AIC", "BIC", "HQ", "max_root"))
  expect_equal(g$AIC, AIC(fit))
  expect_equal(g$max_root, max(stability(fit)))
  a <- augment(fit)
  expect_equal(nrow(a), nrow(yraw))
  expect_true(all(c(".fitted_infl", ".resid_rate") %in% names(a)))
  expect_true(all(is.na(a$.fitted_infl[1:2])))
  expect_equal(
    a$.fitted_infl[-(1:2)] + a$.resid_infl[-(1:2)],
    yraw$infl[-(1:2)]
  )
  expect_error(augment(fit, data = yraw[1:5, ]), "same number of rows")
})

test_that("tidy irf/fevd carry estimates and nested bands", {
  r <- irf(fit, horizon = 4, ci = "mc", nboot = 40)
  td <- tidy(r)
  expect_named(
    td,
    c(
      "response",
      "shock",
      "horizon",
      "estimate",
      "level",
      "conf.low",
      "conf.high"
    )
  )
  expect_s3_class(td$response, "factor")
  expect_equal(nrow(td), 9 * 5 * 2)
  expect_true(all(td$conf.low <= td$conf.high))
  wide <- td[td$level == 0.9, ]
  narrow <- td[td$level == 0.68, ]
  expect_true(all(wide$conf.low <= narrow$conf.low))
  i <- td$response == "rate" &
    td$shock == "infl" &
    td$horizon == 3 &
    td$level == 0.9
  expect_equal(td$estimate[i], r$estimate["rate", "infl", 4])
  expect_equal(tibble::as_tibble(r), td)
  expect_named(
    tidy(r, bands = FALSE),
    c("response", "shock", "horizon", "estimate")
  )
  fe <- tidy(fevd(fit, horizon = 3, ci = "none"))
  expect_equal(nrow(fe), 9 * 4)
})

test_that("tidy hd and forecast", {
  h <- hd(fit)
  td <- tidy(h)
  expect_named(td, c("time", "response", "shock", "contribution", "actual"))
  expect_setequal(levels(td$shock), c(names(yraw), "baseline"))
  tot <- tapply(td$contribution, list(td$time, td$response), sum)
  expect_equal(unname(tot), unname(h$actual), tolerance = 1e-8)
  fc <- forecast(fit, h = 4, nboot = 60)
  tf <- tidy(fc)
  expect_named(
    tf,
    c("horizon", "variable", "estimate", "level", "conf.low", "conf.high")
  )
  expect_equal(nrow(tf), 4 * 3 * 2)
  expect_named(
    tidy(forecast(fit, h = 4, ci = "none")),
    c("horizon", "variable", "estimate")
  )
})
