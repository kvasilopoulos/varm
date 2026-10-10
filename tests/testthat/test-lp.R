n <- nrow(Y3)
p <- 3
H <- 6
lead <- function(x, h) if (h == 0) x else c(x[-(1:h)], rep(NA, h))
lagm <- function(x, l) c(rep(NA, l), x[1:(n - l)])
dfL <- data.frame(Y3)
for (l in 1:p) {
  for (v in colnames(Y3)) {
    dfL[[paste0(v, ".l", l)]] <- lagm(Y3[, v], l)
  }
}
ctr <- paste(
  c("infl", "unemp", setdiff(names(dfL), colnames(Y3))),
  collapse = "+"
)
manual <- function(h, v, z = NULL) {
  d <- dfL
  d$yy <- lead(Y3[, v], h)
  d <- d[(p + 1):(n - h), ]
  suppressWarnings(lm(as.formula(paste("yy ~ rate +", ctr)), d))
}

test_that("LP coefficients and HAC/HC standard errors agree with lm + sandwich", {
  skip_if_not_installed("sandwich")
  f <- lp(
    Y3,
    shock = "rate",
    horizon = H,
    lags = p,
    contemp = c("infl", "unemp")
  )
  f2 <- lp(
    Y3,
    shock = "rate",
    horizon = H,
    lags = p,
    contemp = c("infl", "unemp"),
    se = "hc"
  )
  suppressWarnings(
    for (h in c(0, 3, H)) {
      for (v in colnames(Y3)) {
        fm <- manual(h, v)
        expect_equal(
          unname(f$estimate[h + 1, v]),
          unname(coef(fm)["rate"]),
          tolerance = 1e-8
        )
        expect_equal(
          unname(f$se[h + 1, v]),
          unname(sqrt(sandwich::NeweyWest(
            fm,
            lag = h + 1,
            prewhite = FALSE,
            adjust = FALSE
          )["rate", "rate"])),
          tolerance = 1e-8
        )
        expect_equal(
          unname(f2$se[h + 1, v]),
          unname(sqrt(sandwich::vcovHC(fm, "HC1")["rate", "rate"])),
          tolerance = 1e-8
        )
      }
    }
  )
})

test_that("LP-IV equals manual 2SLS", {
  set.seed(2)
  z <- Y3[, "rate"] + rnorm(n, sd = 0.8)
  z2 <- 0.5 * Y3[, "rate"] + rnorm(n)
  fi <- lp(
    Y3,
    shock = "rate",
    horizon = 4,
    lags = p,
    contemp = c("infl", "unemp"),
    iv = cbind(z, z2),
    se = "hc"
  )
  h <- 2
  d <- dfL
  d$z <- z
  d$z2 <- z2
  d$yy <- lead(Y3[, "infl"], h)
  d <- d[(p + 1):(n - h), ]
  d$rhat <- fitted(lm(as.formula(paste("rate ~ z + z2 +", ctr)), d))
  ss <- lm(as.formula(paste("yy ~ rhat +", ctr)), d)
  expect_equal(
    unname(fi$estimate[h + 1, "infl"]),
    unname(coef(ss)["rhat"]),
    tolerance = 1e-8
  )
})

test_that("LP options and edge cases", {
  f0 <- lp(Y3, "rate", horizon = 0, lags = 1)
  expect_equal(dim(f0$estimate), c(1, 3))
  fa <- lp(Y3, "rate", horizon = 4, lags = 2, lag_augment = TRUE)
  expect_equal(fa$se_type, "hc")
  expect_equal(fa$lags, 3)
  fr <- lp(
    Y3,
    "rate",
    horizon = 4,
    lags = 2,
    responses = c("infl", "unemp"),
    cumulative = "infl"
  )
  expect_equal(colnames(fr$estimate), c("infl", "unemp"))
  fe <- lp(
    Y3,
    rnorm(n),
    horizon = 3,
    lags = 2,
    controls = cbind(q = rnorm(n)),
    dof_adjust = TRUE,
    nw_lag = 2
  )
  expect_equal(dim(fe$bands), c(4, 3, 4))
  expect_error(lp(Y3, "zz"), "zz")
  expect_error(lp(Y3, rnorm(5)), "length")
  expect_error(lp(Y3, "rate", lags = 0), "lags")
  expect_equal(nrow(as.data.frame(fr, conf.int = FALSE)), 5 * 2)
})

test_that("LP and Cholesky VAR agree qualitatively (unit-impact normalisation)", {
  f <- lp(Y3, "rate", horizon = 8, lags = p, contemp = c("infl", "unemp"))
  v <- irf(varm(Y3, p = p), id = id_chol(), horizon = 8, ci = "none")$estimate
  vv <- v["unemp", "rate", ] / v["rate", "rate", 1]
  expect_gt(cor(f$estimate[, "unemp"], vv), 0.8)
})

test_that("variable arguments accept tidyselect and ignore non-numeric columns", {
  ref <- lp(Y3, "rate", horizon = 4, lags = 2, contemp = c("infl", "unemp"))
  df <- tibble::as_tibble(yraw)
  df$date <- as.Date("2000-01-01") + seq_len(nrow(df))
  f <- lp(df, rate, horizon = 4, lags = 2, contemp = c(infl, unemp))
  expect_equal(f$estimate, ref$estimate)
  f2 <- lp(
    df,
    rate,
    horizon = 4,
    lags = 2,
    contemp = !rate,
    cumulative = starts_with("in")
  )
  expect_equal(f2$responses, c("infl", "unemp", "rate"))
  f3 <- lp(df, 3, horizon = 2, lags = 1, responses = all_of(c("rate", "infl")))
  expect_equal(f3$responses, c("rate", "infl"))
  s <- "rate"
  expect_equal(lp(df, !!s, horizon = 2, lags = 1)$shock, "rate")
})

test_that("tidy.lpm and glance.lpm have the documented shape", {
  f <- lp(Y3, "rate", horizon = H, lags = 2, contemp = c("infl", "unemp"))
  td <- tidy(f)
  expect_s3_class(td, "tbl_df")
  expect_named(
    td,
    c(
      "horizon",
      "response",
      "estimate",
      "std.error",
      "level",
      "conf.low",
      "conf.high"
    )
  )
  expect_equal(nrow(td), (H + 1) * 3 * 2)
  expect_s3_class(td$response, "factor")
  expect_setequal(td$level, c(0.68, 0.9))
  expect_true(all(
    td$conf.low <= td$estimate + 1e-12 & td$estimate <= td$conf.high + 1e-12
  ))
  expect_named(
    tidy(f, conf.int = FALSE),
    c("horizon", "response", "estimate", "std.error")
  )
  expect_equal(tibble::as_tibble(f), td)
  gl <- glance(f)
  expect_equal(nrow(gl), 1L)
  expect_true(all(
    c("shock", "horizon", "lags", "se_type", "nobs") %in% names(gl)
  ))
})

test_that("autoplot.lpm returns a ggplot that builds", {
  f <- lp(Y3, "rate", horizon = 4, lags = 2)
  p <- autoplot(f)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_s3_class(plot(f), "ggplot")
})
