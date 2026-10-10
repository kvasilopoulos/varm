# Shared helpers: load varm from source, load BEAR tutorial data, BEAR algorithm replicas in R.
# Usage: source("validation/bear/common.R") from the varm package root.
.libPaths(c(.libPaths(), "/opt/homebrew/lib/R/4.6/site-library"))  # readxl lives in the site library, not the rv one
# Use the installed build in the rv library by default (load_all would recompile src/ and regenerate RcppExports.cpp in a shared working tree);
# set VARM_SRC=1 to test the working-tree sources instead.
if (nzchar(Sys.getenv("VARM_SRC"))) suppressMessages(pkgload::load_all(".", quiet = TRUE)) else suppressMessages(library(varm, lib.loc = "rv/library/4.6/arm64"))
BEAR <- "/Users/VASILK04/Documents/research/external/BEAR-toolbox"
rep_ <- function(name, ok, ...) cat(sprintf("[%s] %s  %s\n", if (is.na(ok)) "NOT-TESTED" else if (ok) "PASS" else "FAIL", name, paste(..., collapse = " ")))

# exampleData.csv (BEARX-tutorials, t1_baseEstimators.m): DOM_GDP, DOM_CPI, STN endogenous, Oil exogenous.
# Estimation span 1975Q1-2014Q4 plus p presample rows (BEAR "longY").
load_bear_data <- function(p = 4, exo = TRUE) {
  d <- read.csv(file.path(BEAR, "BEARX-tutorials/data/exampleData.csv"), check.names = FALSE)
  names(d)[1] <- "Time"
  d <- d[d$Time >= "1974-Q1" & d$Time <= "2014-Q4", ]
  # BEAR longY starts p periods before 1975Q1
  i0 <- which(d$Time == "1975-Q1"); d <- d[(i0 - p):nrow(d), ]
  list(Y = as.matrix(d[, c("DOM_GDP", "DOM_CPI", "STN")]), X = as.matrix(d[, "Oil", drop = FALSE]), time = d$Time)
}

# ---- BEAR olsvar.m ----
bear_olsvar <- function(Y, Xex, const, p) {
  n <- ncol(Y); Tn <- nrow(Y)
  ex <- cbind(if (const) rep(1, Tn), Xex)
  m <- if (is.null(ex)) 0 else ncol(ex)
  if (m) ex <- ex[(p + 1):Tn, , drop = FALSE]
  lag <- do.call(cbind, lapply(1:p, function(l) Y[(p + 1 - l):(Tn - l), , drop = FALSE]))
  X <- cbind(lag, ex); Yd <- Y[(p + 1):Tn, , drop = FALSE]
  B <- solve(crossprod(X), crossprod(X, Yd)); E <- Yd - X %*% B
  k <- ncol(X)
  list(B = B, Sigma = crossprod(E) / (nrow(X) - k), X = X, Y = Yd, E = E, n = n, m = m, k = k, T = nrow(X))
}
# ---- BEAR arloop.m: AR(p)+const residual variance, divisor T-k with k = p + 1 ----
bear_arvar <- function(Y, const, p) vapply(seq_len(ncol(Y)), function(i) bear_olsvar(Y[, i, drop = FALSE], NULL, const, p)$Sigma[1, 1], 0)

# ---- BEAR nwprior.m (prior 21 'ar') + dummies + nwpost.m ----
bear_nw <- function(Y, Xex, p, const = TRUE, lambda1 = .1, lambda3 = 1, lambda4 = 100, ar = .8,
                    sumcoef = NULL, initobs = NULL, sigma = c("ar", "eye"), ini0 = 0) {
  # ini0 = 0: BEAR (dummy mean over the p PRE-sample rows); ini0 = p mimics varm (first p rows of the regression sample)
  sigma <- match.arg(sigma)
  o <- bear_olsvar(Y, Xex, const, p); n <- o$n; m <- o$m; k <- o$k
  arv <- bear_arvar(Y, const, p)
  B0 <- matrix(0, k, n); B0[cbind(1:n, 1:n)] <- ar
  phi0 <- numeric(k)
  for (j in 1:p) for (i in 1:n) phi0[(j - 1) * n + i] <- (1 / arv[i]) * (lambda1 / j^lambda3)^2
  if (m) phi0[(k - m + 1):k] <- (lambda1 * lambda4)^2
  alpha0 <- n + 2
  S0 <- if (sigma == "ar") (alpha0 - n - 1) * diag(arv) else diag(n)
  Yd <- o$Y; Xd <- o$X; ini <- Y[ini0 + 1:p, , drop = FALSE]
  if (!is.null(sumcoef)) {  # dummies.SumCoeff
    dy <- diag(colMeans(ini) / sumcoef); Yd <- rbind(Yd, dy)
    Xd <- rbind(Xd, cbind(do.call(cbind, rep(list(dy), p)), matrix(0, n, m)))
  }
  if (!is.null(initobs)) {  # dummies.InitialObs
    dy <- matrix(colMeans(ini) / initobs, 1)
    ex0 <- if (m) colMeans(cbind(if (const) rep(1, p), if (!is.null(Xex)) Xex[ini0 + 1:p, , drop = FALSE])) / initobs else NULL
    Yd <- rbind(Yd, dy); Xd <- rbind(Xd, c(do.call(cbind, rep(list(dy), p)), ex0))
  }
  Tn <- nrow(Yd); ip <- diag(1 / phi0, k); ipb <- ip + crossprod(Xd)
  phibar <- solve(ipb); Bbar <- phibar %*% (ip %*% B0 + crossprod(Xd, Yd))
  Sbar <- crossprod(Yd) + S0 + t(B0) %*% ip %*% B0 - t(Bbar) %*% ipb %*% Bbar
  Sbar <- (Sbar + t(Sbar)) / 2
  list(Bbar = Bbar, phibar = phibar, Sbar = Sbar, alphabar = Tn + alpha0, arvar = arv, phi0 = phi0, B0 = B0, S0 = S0,
       ESigma = Sbar / (Tn + alpha0 - n - 1))
}

# ---- BEAR irfsim.m (recursion) : returns n x n x H array [resp, shock, h] of reduced-form IRFs and orthogonalised ----
bear_irfsim <- function(B, D, n, m, p, H) {
  ir <- array(0, c(n, n, H))
  for (i in 1:n) {
    Y <- matrix(0, p, n); Y[p, i] <- 1
    for (j in 1:(H - 1)) {
      x <- c(t(Y[nrow(Y):(nrow(Y) - p + 1), , drop = FALSE])); x <- c(x, rep(0, m))
      Y <- rbind(Y, x %*% B)
    }
    ir[, i, ] <- t(Y[p:nrow(Y), ])
  }
  ort <- array(0, c(n, n, H)); for (h in 1:H) ort[, , h] <- ir[, , h] %*% D
  list(irf = ir, ort = ort)
}

# ---- read BEAR-shipped expected results (tests/bear5/*.mat) via the pure-python reader ----
bear_mat <- function(file, vars) {
  out <- tempfile(fileext = ".json")
  st <- system2("python3", c(shQuote("validation/bear/matread.py"), shQuote(file.path(BEAR, "tests/bear5", file)), "--json", out, vars))
  if (st != 0) stop("matread failed")
  j <- jsonlite::fromJSON(out, simplifyVector = FALSE)
  conv <- function(n) {
    if (!is.null(n$cell)) { r <- lapply(n$cell, conv); dim(r) <- unlist(n$dim); r }
    else if (!is.null(n$x)) { x <- as.numeric(unlist(lapply(n$x, function(z) if (is.null(z)) NA else z))); if (length(n$dim) == 2) matrix(x, n$dim[[1]], n$dim[[2]]) else x }
    else if (!is.null(n$str)) n$str else NULL
  }
  lapply(j, conv)
}

# ---- BEAR mprior.m / mpost.m (Minnesota; prior 11/12/13) ----
bear_minnesota <- function(Y, p, const = TRUE, ar = .8, l1 = .1, l2 = .5, l3 = 1, l4 = 100, prior = 12) {
  o <- bear_olsvar(Y, NULL, const, p); n <- o$n; m <- o$m; k <- o$k; q <- n * k
  arv <- bear_arvar(Y, const, p)
  beta0 <- numeric(q); for (i in 1:n) beta0[(i - 1) * k + i] <- ar
  om <- numeric(q)
  for (i in 1:n) for (j in 1:p) for (kk in 1:n) {
    ix <- (i - 1) * k + (j - 1) * n + kk
    om[ix] <- if (kk == i) (l1 / j^l3)^2 else (arv[i] / arv[kk]) * ((l1 * l2) / j^l3)^2
  }
  for (i in 1:n) for (jj in seq_len(m)) om[i * k - m + jj] <- arv[i] * (l1 * l4)^2
  sigma <- switch(as.character(prior), "11" = diag(arv), "12" = diag(diag(o$Sigma)), "13" = o$Sigma)
  iom <- diag(1 / om); isg <- solve(sigma)
  iob <- iom + kronecker(isg, crossprod(o$X)); omb <- solve(iob)
  bb <- omb %*% (iom %*% beta0 + kronecker(isg, t(o$X)) %*% c(o$Y))
  list(beta0 = beta0, omega0 = diag(om), sigma = sigma, betabar = bb, omegabar = omb, arvar = arv)
}
