# BEAR vs varm: OLS VAR, AR variances, Normal-Wishart posterior (incl. dummies), run from package root.
source("validation/bear/common.R")
set.seed(1)
p <- 4; d <- load_bear_data(p); Y <- d$Y; X <- d$X; cat("T(longY)=", nrow(Y), "\n")

## T1 OLS VAR (const + exogenous Oil)
o <- bear_olsvar(Y, X, TRUE, p); f <- varm(Y, p, const = TRUE, exog = X)
rn <- rownames(f$B); idx <- c(1:(3 * p), which(rn == "const"), which(rn == "Oil")) # BEAR order: lags, const, Oil
rep_("T1 OLS B (const+exo)", max(abs(o$B - f$B[idx, ])) < 1e-10, "maxdiff", format(max(abs(o$B - f$B[idx, ])), digits = 3))
rep_("T1 OLS Sigma (T-k divisor)", max(abs(o$Sigma - f$Sigma)) < 1e-10, "maxdiff", format(max(abs(o$Sigma - f$Sigma)), digits = 3))

## T2 AR residual variances (Minnesota/NW scale)
a_b <- bear_arvar(Y, TRUE, p); a_v <- varm:::ar_resid_var(Y, p)
cat("  BEAR arvar", format(a_b, digits = 5), "\n  varm arvar", format(a_v, digits = 5), "\n  ratio", format(a_v / a_b, digits = 5),
    "  expected (T-k)/(n-p) =", (nrow(Y) - p - (p + 1)) / (nrow(Y) - p), "\n")
rep_("T2 AR variance identical", max(abs(a_b / a_v - 1)) < 1e-8, "(convention: BEAR /(T-p-1), varm /(n-p))")

## T3/T4 NW posterior, hyperparameters mapped: lambda=lambda1, alpha=lambda3, psi=BEAR arvar, delta=ar, var=(l1*l4)^2
chk <- function(nm, Y, X, sc = NULL, io = NULL, n = 40000, ex = TRUE, ini0 = 0) {
  nw <- bear_nw(Y, if (ex) X else NULL, p, sumcoef = sc, initobs = io, ini0 = ini0)
  nm <- paste(nm, if (ini0) "[BEAR dummies built from regression-sample rows, i.e. varm convention]" else "")
  pr <- prior_minnesota(lambda = .1, alpha = 2, soc = if (is.null(sc)) NULL else c(sc, 1, 1e-4, 50),
                        sur = if (is.null(io)) NULL else c(io, 1, 1e-4, 50), psi = nw$arvar, delta = .8,
                        hyper = character(0), var = (.1 * 100)^2)
  fit <- bvarm(Y, p, prior = pr, ndraw = n, exog = if (ex) X else NULL)
  B <- fit$draws$B; S <- fit$draws$Sigma; rn <- dimnames(B)[[1]]
  ix <- if (ex) c(1:(3 * p), which(rn == "const"), which(rn == "Oil")) else seq_len(dim(B)[1])
  Bm <- apply(B, c(1, 2), mean)[ix, ]; Bs <- apply(B, c(1, 2), sd)[ix, ] / sqrt(n)
  zB <- max(abs((Bm - nw$Bbar) / Bs)); relB <- max(abs(Bm - nw$Bbar)) 
  Sm <- apply(S, c(1, 2), mean); zS <- max(abs((Sm - nw$ESigma) / (apply(S, c(1, 2), sd) / sqrt(n))))
  # analytic posterior sd of B coefficient: sqrt(phibar_ii * Sbar_jj/(alpha-K-1))
  sdA <- sqrt(outer(diag(nw$phibar), diag(nw$Sbar) / (nw$alphabar - 3 - 1))); sdD <- apply(B, c(1, 2), sd)[ix, ]
  rep_(paste(nm, "E[B] z"), zB < 4.5, "max|z|=", round(zB, 2), " max abs diff", format(relB, digits = 3))
  rep_(paste(nm, "E[Sigma] z"), zS < 4.5, "max|z|=", round(zS, 2))
  rep_(paste(nm, "sd(B) ratio"), max(abs(sdD / sdA - 1)) < 0.03, "max|ratio-1|=", round(max(abs(sdD / sdA - 1)), 4))
  invisible(list(fit = fit, nw = nw))
}
chk("T3 NW no-dummy, const+Oil", Y, X)
chk("T4a NW + SumCoeff(0.1), const+Oil", Y, X, sc = .1)
chk("T4b NW + InitialObs(1e-3) no exog", Y, X, io = 1e-3, ex = FALSE)
chk("T4c NW + InitialObs(1e-3) with exog Oil", Y, X, io = 1e-3, ex = TRUE)

cat("--- diagnosis: replicate BEAR dummies with the mean of the first p rows of the *regression* sample (what varm does)\n")

## T11 marginal likelihood: BEAR nwmlik.m / nwmlikdata.m (log ML of the NW prior, used by BEAR hyperparameter optimisation) vs varm$logml
mgamma <- function(a, n) { j <- 1:n; n * (n - 1) / 4 * log(pi) + sum(lgamma(a + (1 - j) / 2)) }
bear_nwmlik <- function(Y, Xex, p, sc = NULL, io = NULL, ini0 = 0, dummy_update = TRUE, ...) {
  nw <- bear_nw(Y, Xex, p, sumcoef = sc, initobs = io, ini0 = ini0, ...)
  o <- bear_olsvar(Y, Xex, TRUE, p); n <- 3; k <- o$k; Tn <- o$T
  data_ml <- function(X, T, Sbar, abar, phi0, S0, a0) {
    t1 <- -n * T / 2 * log(pi) + mgamma(abar / 2, n) - mgamma(a0 / 2, n); t2 <- -T / 2 * sum(log(diag(S0)))
    ev <- eigen(diag(sqrt(phi0)) %*% crossprod(X) %*% diag(sqrt(phi0)), symmetric = TRUE)$values; t3 <- -n / 2 * sum(log(1 + ev))
    Fs <- diag(1 / sqrt(diag(S0))); ev2 <- eigen(Fs %*% (Sbar - S0) %*% Fs, symmetric = TRUE)$values; t4 <- -abar / 2 * sum(log(1 + ev2))
    t1 + t2 + t3 + t4 }
  a0 <- n + 2; tot <- nrow(o$X) + (if (!is.null(sc)) n else 0) + (if (!is.null(io)) 1 else 0)
  # rebuild the dummy-augmented X as bear_nw does
  ini <- Y[ini0 + 1:p, , drop = FALSE]; Xd <- o$X; Yd <- o$Y; m <- o$m
  Xdum <- matrix(0, 0, k); Ydum <- matrix(0, 0, n)
  if (!is.null(sc)) { dy <- diag(colMeans(ini) / sc); Ydum <- rbind(Ydum, dy); Xdum <- rbind(Xdum, cbind(do.call(cbind, rep(list(dy), p)), matrix(0, n, m))) }
  if (!is.null(io)) { dy <- matrix(colMeans(ini) / io, 1); ex0 <- colMeans(cbind(rep(1, p), if (!is.null(Xex)) Xex[ini0 + 1:p, , drop = FALSE])) / io
    Ydum <- rbind(Ydum, dy); Xdum <- rbind(Xdum, c(do.call(cbind, rep(list(dy), p)), ex0)) }
  Xa <- rbind(o$X, Xdum); Ya <- rbind(o$Y, Ydum)
  lm_all <- data_ml(Xa, nrow(Xa), nw$Sbar, nw$alphabar, nw$phi0, nw$S0, a0)
  if (nrow(Xdum)) { ip <- diag(1 / nw$phi0); ipb <- ip + crossprod(Xdum); pb <- solve(ipb); Bb <- pb %*% (ip %*% nw$B0 + crossprod(Xdum, Ydum))
    Sb <- crossprod(Ydum) + nw$S0 + t(nw$B0) %*% ip %*% nw$B0 - t(Bb) %*% ipb %*% Bb; Sb <- (Sb + t(Sb)) / 2
    # dummy_update = FALSE: residuals of the dummies around the PRIOR mean B0, no posterior update (what varm's C++ ev(..., est=false) does)
    if (!dummy_update) { R <- Ydum - Xdum %*% nw$B0; Sb <- nw$S0 + crossprod(R) }
    lm_all <- lm_all - data_ml(Xdum, nrow(Xdum), Sb, nrow(Xdum) + a0, nw$phi0, nw$S0, a0) }
  list(logml = lm_all, nw = nw)
}
chk_ml <- function(nm, sc = NULL, io = NULL, ex = TRUE) {
  Xe <- if (ex) X else NULL
  b <- bear_nwmlik(Y, Xe, p, sc, io, ini0 = 0)   # dummies from the p presample obs (BEAR convention, now also varm)
  pr <- prior_minnesota(lambda = .1, alpha = 2, soc = if (is.null(sc)) NULL else c(sc, 1, 1e-4, 50), sur = if (is.null(io)) NULL else c(io, 1, 1e-4, 50),
                        psi = b$nw$arvar, delta = .8, hyper = character(0), var = 100)
  fit <- bvarm(Y, p, prior = pr, ndraw = 10, exog = Xe)
  rep_(paste("T11", nm, "log ML varm vs BEAR nwmlik"), abs(fit$logml - b$logml) < 1e-6, "varm", format(fit$logml, digits = 10), "BEAR-formula", format(b$logml, digits = 10))
  if (!is.null(sc) || !is.null(io)) { b2 <- bear_nwmlik(Y, Xe, p, sc, io, ini0 = 0, dummy_update = FALSE)
    rep_(paste("T11'", nm, "varm == BEAR formula with dummy block evaluated at prior mean (diagnosis)"), abs(fit$logml - b2$logml) < 1e-6, "diff", format(fit$logml - b2$logml, digits = 3)) }
}
chk_ml("NW no dummies"); chk_ml("NW + SumCoeff(0.1)", sc = .1); chk_ml("NW + SumCoeff + InitialObs (no exog)", sc = .1, io = .001, ex = FALSE)

## T11'' independent check of which dummy-ML convention is right: p(Y | dummies) computed directly as the NW marginal likelihood of Y under the
## prior obtained by updating the base NW prior with the dummy observations (general NW formula, no diagonal shortcuts).
nw_logml_general <- function(Y, X, B0, Phi0, S0, a0) {
  n <- ncol(Y); T <- nrow(Y); Pi <- solve(Phi0); ipb <- Pi + crossprod(X); Phib <- solve(ipb); Bb <- Phib %*% (Pi %*% B0 + crossprod(X, Y))
  Sb <- crossprod(Y) + S0 + t(B0) %*% Pi %*% B0 - t(Bb) %*% ipb %*% Bb; Sb <- (Sb + t(Sb)) / 2; ab <- a0 + T
  -n * T / 2 * log(pi) + mgamma(ab / 2, n) - mgamma(a0 / 2, n) + a0 / 2 * determinant(S0)$modulus[1] - ab / 2 * determinant(Sb)$modulus[1] +
    n / 2 * (determinant(Phib)$modulus[1] - determinant(Phi0)$modulus[1])
}
{
  sc <- .1; nw <- bear_nw(Y, X, p, sumcoef = sc, ini0 = 0); o <- bear_olsvar(Y, X, TRUE, p); n <- 3; m <- o$m; k <- o$k
  ini <- Y[p + 1:p, , drop = FALSE]; dy <- diag(colMeans(ini) / sc); Xd <- cbind(do.call(cbind, rep(list(dy), p)), matrix(0, n, m))
  Pi0 <- diag(1 / nw$phi0); ipd <- Pi0 + crossprod(Xd); Phid <- solve(ipd); Bd <- Phid %*% (Pi0 %*% nw$B0 + crossprod(Xd, dy))
  Sd <- crossprod(dy) + nw$S0 + t(nw$B0) %*% Pi0 %*% nw$B0 - t(Bd) %*% ipd %*% Bd; Sd <- (Sd + t(Sd)) / 2
  direct <- nw_logml_general(o$Y, o$X, Bd, Phid, Sd, n + 2 + n)
  full <- nw_logml_general(rbind(o$Y, dy), rbind(o$X, Xd), nw$B0, diag(nw$phi0), nw$S0, n + 2) - nw_logml_general(dy, Xd, nw$B0, diag(nw$phi0), nw$S0, n + 2)
  rep_("T11'' log p(Y|SumCoeff dummies): direct sequential-prior formula == BEAR difference-of-MLs", abs(direct - full) < 1e-6, "direct", format(direct, digits = 10), "BEAR-type", format(full, digits = 10))
  cat("   (varm reported -595.5317 for the same quantity; correct value is", format(direct, digits = 10), ")\n")
}
