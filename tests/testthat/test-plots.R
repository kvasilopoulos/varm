fit <- varm(yraw, p = 2)
grDevices::pdf(NULL)

smoke <- function(p) {
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_true(length(b$data) >= 1)
  invisible(b)
}

test_that("theme and palette", {
  expect_s3_class(theme_varm(), "theme")
  expect_length(varm_palette(3), 3)
  expect_match(varm_palette(8), "^#")
  expect_warning(varm_palette(9), "recycling")
})

test_that("irf / fevd autoplot", {
  r <- irf(fit, horizon = 5, ci = "mc", nboot = 40)
  b <- smoke(ggplot2::autoplot(r))
  expect_length(unique(b$layout$layout$PANEL), 9)
  b2 <- smoke(ggplot2::autoplot(r, shocks = rate, responses = c(infl, unemp)))
  expect_equal(nrow(b2$layout$layout), 2)
  smoke(ggplot2::autoplot(irf(fit, horizon = 5, ci = "none")))
  smoke(ggplot2::autoplot(fevd(fit, horizon = 5, ci = "mc", nboot = 30)))
  expect_error(ggplot2::autoplot(r, shocks = zzz), "unknown shock")
})

test_that("forecast autoplot: history and fan, with and without dates", {
  fc <- forecast(fit, h = 6, nboot = 60)
  smoke(ggplot2::autoplot(fc))
  smoke(ggplot2::autoplot(fc, variables = rate, n_history = 10))
  smoke(ggplot2::autoplot(forecast(fit, h = 3, ci = "none")))
  d <- cbind(
    date = seq(as.Date("1990-01-01"), by = "quarter", length.out = nrow(yraw)),
    yraw
  )
  smoke(ggplot2::autoplot(forecast(varm(d, p = 2), h = 4, nboot = 30)))
})

test_that("hd autoplot", {
  smoke(ggplot2::autoplot(hd(fit)))
  smoke(ggplot2::autoplot(hd(fit), responses = c(rate, infl)))
})

test_that("draws autoplot: trace and density", {
  d <- draws(fit, "mc", 50)
  smoke(ggplot2::autoplot(d))
  smoke(ggplot2::autoplot(d, type = "density", terms = c(rate.l1, const)))
  expect_error(ggplot2::autoplot(d, type = "x"))
})

test_that("plot() wrappers return the ggplot invisibly", {
  p <- withVisible(plot(irf(fit, horizon = 3, ci = "none")))
  expect_false(p$visible)
  expect_s3_class(p$value, "ggplot")
  expect_s3_class(plot(forecast(fit, h = 3, nboot = 20)), "ggplot")
  expect_s3_class(plot(hd(fit)), "ggplot")
  expect_s3_class(plot(draws(fit, "mc", 20)), "ggplot")
})
