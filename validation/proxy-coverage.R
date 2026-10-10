# Monte Carlo coverage of proxy-SVAR 90% intervals: wild bootstrap vs Jentsch-Lunsford
# moving block bootstrap vs Montiel Olea-Stock-Watson (MOSW) robust sets.
# Not run in tests. Usage: Rscript validation/proxy-coverage.R [reps] [nboot] [T]
# DGP: bivariate VAR(1), iid N(0,1) structural shocks, proxy m = lambda * e1 + sd_v * v.
args <- as.numeric(commandArgs(TRUE))
reps <- if (length(args) >= 1) args[1] else 200
nboot <- if (length(args) >= 2) args[2] else 199
Tn <- if (length(args) >= 3) args[3] else 200
pkgload::load_all(".", quiet = TRUE)

A <- matrix(c(0.5, 0.2, 0.1, 0.4), 2, 2)
B <- matrix(c(1, 0.5, 0.3, 1), 2, 2)
H <- 8L
hs <- c(1L, 4L, 8L) # horizons reported (impact of y2 is also reported, h = 0)
Th <- varm:::cpp_ma(t(A), 2L, 1L, H)
truth <- t(vapply(0:H, function(h) drop(Th[,, h + 1] %*% B[, 1] / B[1, 1]), numeric(2)))

designs <- list(
  strong = c(lambda = 1, sd_v = 0.5),
  weak = c(lambda = 0.1, sd_v = 1)
)

one_rep <- function(seed, des) {
  set.seed(seed)
  e <- matrix(rnorm((Tn + 50) * 2), ncol = 2)
  u <- e %*% t(B)
  y <- matrix(0, Tn + 50, 2)
  for (t in 2:(Tn + 50)) y[t, ] <- A %*% y[t - 1, ] + u[t, ]
  keep <- 51:(Tn + 50)
  y <- y[keep, ]
  m <- des[["lambda"]] * e[keep, 1] + des[["sd_v"]] * rnorm(Tn)
  fit <- varm(y, p = 1)
  res <- list()
  cov_boot <- function(r) {
    lo <- r$bands[, 2, "0.05"]
    hi <- r$bands[, 2, "0.95"]
    cbind(cover = lo <= truth[, 2] & truth[, 2] <= hi, width = hi - lo)
  }
  res$wild <- cov_boot(varm_proxy(fit, m, horizon = H, nboot = nboot, level = 0.9, method = "wild"))
  res$mbb <- cov_boot(varm_proxy(fit, m, horizon = H, nboot = nboot, level = 0.9, method = "mbb"))
  rb <- varm_proxy_robust(fit, m, horizon = H, level = 0.9)
  cs <- rb$case[, 2, 1]
  lo <- rb$lower[, 2, 1]
  hi <- rb$upper[, 2, 1]
  gap <- rb$gap[, 2, 1, ]
  cov <- ifelse(
    cs == 1, lo <= truth[, 2] & truth[, 2] <= hi,
    ifelse(cs == 2, !(gap[, 1] < truth[, 2] & truth[, 2] < gap[, 2]),
      ifelse(cs == 4, TRUE, FALSE))
  )
  res$mosw <- cbind(cover = cov, width = ifelse(cs == 1, hi - lo, Inf))
  res$F <- rb$fstat[["F_hac"]]
  res
}

run <- function(des_name) {
  des <- designs[[des_name]]
  out <- parallel::mclapply(seq_len(reps), function(i) one_rep(1000 + i, des), mc.cores = max(1, parallel::detectCores() - 1))
  ok <- !vapply(out, inherits, NA, "try-error")
  out <- out[ok]
  tab <- do.call(rbind, lapply(c("wild", "mbb", "mosw"), function(mth) {
    cv <- sapply(out, function(o) o[[mth]][, "cover"])
    wd <- sapply(out, function(o) o[[mth]][, "width"])
    data.frame(
      design = des_name, method = mth, reps = length(out),
      h = c(0L, hs), cover = rowMeans(cv)[c(1, hs + 1)],
      median_width = apply(wd, 1, median)[c(1, hs + 1)]
    )
  }))
  cat("\n== ", des_name, ": median HAC first-stage F = ",
      round(median(sapply(out, `[[`, "F")), 1), " ==\n", sep = "")
  print(tab, row.names = FALSE, digits = 3)
  invisible(tab)
}

cat("T =", Tn, " reps =", reps, " nboot =", nboot, " nominal = 0.90; response of y2 per unit y1 impact\n")
res <- lapply(names(designs), run)
