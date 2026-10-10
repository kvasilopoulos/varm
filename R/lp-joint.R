# Cross-horizon covariance helpers -------------------------------------------

# Newey-West long-run variance of the columns of G (n x q), Bartlett kernel,
# truncation L
nw_lrv <- function(G, L) {
  G <- as.matrix(G)
  n <- nrow(G)
  out <- crossprod(G)
  for (l in seq_len(min(L, n - 1L))) {
    C <- crossprod(G[(l + 1):n, , drop = FALSE], G[1:(n - l), , drop = FALSE])
    out <- out + (1 - l / (L + 1)) * (C + t(C))
  }
  out
}

# Wald statistic b' V^-1 b, NA if V is numerically singular; zero-variance
# (normalised) horizons are dropped
wald_stat <- function(b, V) {
  ok <- diag(V) > 1e-12 * max(diag(V))
  b <- b[ok]
  V <- V[ok, ok, drop = FALSE]
  st <- tryCatch(drop(crossprod(b, solve(V, b))), error = function(e) NA_real_)
  p <- if (is.na(st)) {
    NA_real_
  } else {
    stats::pchisq(st, length(b), lower.tail = FALSE)
  }
  c(stat = st, df = length(b), p = p)
}

#' Local projections with joint cross-horizon covariance
#'
#' Same regressions as [lp()] but additionally returns, for every response, the covariance matrix
#' of the whole impulse-response vector \eqn{(\hat\beta_0,\dots,\hat\beta_H)'},
#' \deqn{\widehat{\mathrm{Cov}}(\hat\beta_h,\hat\beta_{h'}) = \sum_{l=-L}^{L} w_l \sum_t
#'   q_{t,h} e_{t+h}\, q_{t-l,h'} e_{t-l+h'},\quad q_{t,h} = x_{t,h}'(X_h'X_h)^{-1}e_1,}
#' with Bartlett weights \eqn{w_l = 1-|l|/(L+1)} and \eqn{L = \max(h,h')+1} (or a common `nw_lag`).
#' The diagonal reproduces the standard errors of [lp()]. The engine is `cpp_lp_cov` (C++).
#'
#' **Economic intuition.** Local-projection estimates at neighbouring horizons are built from
#' overlapping data and highly correlated, so reading a pointwise band as a statement about the
#' whole response path overstates what the data say. The joint covariance is the input for
#' simultaneous inference ([lp_supt()]), joint hypothesis tests (e.g. "the response is zero at all
#' horizons"), and smoothing. *Pitfall:* with the default horizon-specific bandwidth the full matrix is
#' not guaranteed positive semi-definite; use a common `nw_lag` if you need to invert it.
#'
#' @inheritParams lp
#' @return an `lp_joint` object (also an `lpm`) with the extra element `$vcov`: a named list (one per
#'   response) of \eqn{(H+1)\times(H+1)} covariance matrices.
#' @family local projections
#' @examples
#' fit <- lp_joint(yraw, shock = rate, horizon = 8, lags = 2, contemp = c(infl, unemp))
#' sqrt(diag(fit$vcov[[1]])) - fit$se[, 1] # zero
#' @export
lp_joint <- function(
  Y,
  shock,
  horizon = 20L,
  lags = 4L,
  contemp = NULL,
  controls = NULL,
  iv = NULL,
  responses = NULL,
  cumulative = NULL,
  se = c("nw", "hc"),
  nw_lag = NULL,
  lag_augment = FALSE,
  level = c(0.68, 0.90),
  dof_adjust = FALSE
) {
  spec <- lp_spec(
    Y,
    rlang::enquo(shock),
    rlang::enquo(contemp),
    rlang::enquo(responses),
    rlang::enquo(cumulative),
    horizon = horizon,
    lags = lags,
    controls = controls,
    iv = iv,
    se = match.arg(se),
    nw_lag = nw_lag,
    lag_augment = lag_augment,
    level = level,
    dof_adjust = dof_adjust,
    call = match.call()
  )
  fit <- lp_fit(spec)
  r <- do.call(cpp_lp_cov, spec$engine)
  hh <- 0:spec$H
  fit$vcov <- lapply(seq_along(fit$responses), function(k) {
    matrix(r$vcov[,, k], length(hh), dimnames = list(hh, hh))
  })
  names(fit$vcov) <- fit$responses
  class(fit) <- c("lp_joint", class(fit))
  fit
}

# critical value of max_h |z_h|, z ~ N(0, cov2cor(V))
supt_crit <- function(V, level, nsim, seed = NULL) {
  ok <- diag(V) > 1e-14
  R <- diag(nrow(V))
  R[ok, ok] <- stats::cov2cor(V[ok, ok, drop = FALSE])
  ev <- eigen(R, symmetric = TRUE)
  A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(ev$values))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  z <- matrix(stats::rnorm(nsim * nrow(V)), nsim) %*% t(A)
  unname(stats::quantile(apply(abs(z), 1, max), level))
}


#' Sup-t simultaneous confidence bands for impulse responses
#'
#' Montiel Olea & Plagborg-Moller (2019) sup-t bands. With \eqn{\hat\sigma_h} the standard error of
#' \eqn{\hat\beta_h} and \eqn{\hat R} the correlation matrix implied by the joint covariance, the band is
#' \deqn{\hat\beta_h \pm \hat c_{1-\alpha}\hat\sigma_h,\qquad \hat c_{1-\alpha} = (1-\alpha)\text{-quantile of }
#' \max_h |z_h|,\ z\sim N(0,\hat R),}
#' computed by simulation. The band covers the entire response path with probability \eqn{1-\alpha}
#' (asymptotically), is narrower than Bonferroni and is never narrower than the pointwise band.
#'
#' **Economic intuition.** A pointwise 90% band is a statement about one horizon at a time; across 20
#' horizons the chance that the true path leaves it somewhere is far above 10%. Use sup-t bands
#' when you want to claim that a whole shape (a hump, a sign, "no response at any horizon") is
#' supported. Because LP estimates are strongly positively correlated across horizons the sup-t
#' critical value is usually only modestly above 1.64 (typically 2.0-2.5), not the 3+ of
#' Bonferroni. *Pitfall:* it needs a good joint covariance (heavy persistence makes HAC covariances
#' unreliable at long horizons) and horizons with zero variance (normalised impact) are ignored.
#'
#' @param object an object with a `$vcov` element: [lp_joint()], [lp_smooth()] or [lp_panel()].
#' @param level simultaneous coverage level.
#' @param nsim number of simulation draws.
#' @param seed optional RNG seed.
#' @return a tibble with columns `response`, `horizon`, `estimate`, `std.error`,
#'   `pw_low`, `pw_high` (pointwise at `level`), `conf.low`, `conf.high` (sup-t);
#'   the critical values are in `attr(, "crit")`. [tidy()][tidy.lp_supt()] stacks the two
#'   band types and [autoplot()] compares them.
#' @family local projections
#' @examples
#' fit <- lp_joint(yraw, shock = rate, horizon = 8, lags = 2, contemp = c(infl, unemp))
#' b <- lp_supt(fit, level = 0.9, seed = 1)
#' attr(b, "crit")
#' @export
lp_supt <- function(object, level = 0.90, nsim = 10000L, seed = NULL) {
  if (is.null(object$vcov)) {
    rlang::abort(
      "`object` has no joint covariance: use lp_joint(), lp_smooth() or lp_panel()"
    )
  }
  est <- as.matrix(object$estimate)
  V <- object$vcov
  if (is.matrix(V)) {
    V <- list(V)
  }
  zp <- stats::qnorm(1 - (1 - level) / 2)
  one <- function(k) {
    sd_ <- sqrt(pmax(diag(V[[k]]), 0))
    cc <- supt_crit(V[[k]], level, nsim, seed)
    tibble::tibble(
      response = colnames(est)[k],
      horizon = seq_len(nrow(est)) - 1L,
      estimate = est[, k],
      std.error = sd_,
      pw_low = est[, k] - zp * sd_,
      pw_high = est[, k] + zp * sd_,
      conf.low = est[, k] - cc * sd_,
      conf.high = est[, k] + cc * sd_,
      crit = cc
    )
  }
  out <- do.call(rbind, lapply(seq_len(ncol(est)), one))
  crit <- stats::setNames(
    out$crit[match(colnames(est), out$response)],
    colnames(est)
  )
  structure(
    out[setdiff(names(out), "crit")],
    crit = crit,
    level = level,
    class = c("lp_supt", class(out))
  )
}
