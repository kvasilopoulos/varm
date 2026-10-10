#' Direct multi-step forecasts from local projections
#'
#' For each horizon \eqn{h=1,\dots,H} regresses \eqn{y_{t+h}} on \eqn{x_t=(1,y_t',\dots,y_{t-p+1}')'} and
#' forecasts \eqn{\hat y_{T+h}=x_T'\hat B_h}. The forecast standard error is
#' \eqn{\hat\sigma_h\sqrt{1+x_T'(X_h'X_h)^{-1}x_T}} (the direct-regression residual has an MA(h-1)
#' structure but its unconditional variance is what matters for the forecast error).
#'
#' **Economic intuition.** Iterated (VAR) forecasts are efficient if the VAR is correctly
#' specified but compound any misspecification with the horizon; direct forecasts
#' estimate each horizon's conditional mean separately, which is robust to misspecification but noisier
#' (Marcellino, Stock & Watson 2006). Prefer direct forecasts at long horizons with long samples
#' and a possibly nonlinear/mis-specified dynamics, iterated forecasts when the sample is short. At
#' \eqn{h=1} the two coincide exactly. *Pitfall:* forecasts at different horizons are not
#' constrained to be consistent with each other and can be jagged.
#'
#' @param Y data frame / tibble / matrix (T x K); non-numeric columns are dropped.
#' @param horizon maximum forecast horizon.
#' @param lags number of lags \eqn{p}.
#' @param level confidence levels for the bands.
#' @return an `lp_forecast` object: `$mean`, `$se` (`H x K`), `$bands`
#'   (`H x K x 2*length(level)`). [tidy()][tidy.lp_forecast()] returns a tibble.
#' @family local projections
#' @examples
#' lp_forecast(yraw, horizon = 4, lags = 2)
#' @export
lp_forecast <- function(Y, horizon = 8L, lags = 4L, level = c(0.68, 0.90)) {
  Y <- prep_data(Y)$Y
  n <- nrow(Y)
  K <- ncol(Y)
  H <- as.integer(horizon)
  p <- as.integer(lags)
  # regressors at time t: (1, y_t, ..., y_{t-p+1}) = lags 0..p-1
  xT <- c(1, as.vector(t(Y[n:(n - p + 1L), , drop = FALSE])))
  X <- cbind(1, Y, if (p > 1L) lag_matrix(Y, p - 1L))
  dn <- list(seq_len(H), colnames(Y))
  est <- se <- matrix(NA_real_, H, K, dimnames = dn)
  for (h in seq_len(H)) {
    ix <- p:(n - h)
    Xh <- X[ix, , drop = FALSE]
    fit <- stats::lm.fit(Xh, Y[ix + h, , drop = FALSE])
    s2 <- colSums(fit$residuals^2) / (length(ix) - ncol(Xh))
    lev <- drop(crossprod(xT, chol2inv(chol(crossprod(Xh))) %*% xT))
    est[h, ] <- drop(crossprod(xT, fit$coefficients))
    se[h, ] <- sqrt(s2 * (1 + lev))
  }
  probs <- level_probs(level)
  structure(
    list(
      mean = est,
      se = se,
      bands = band_array(est, se, probs),
      probs = probs,
      level = level,
      h = H,
      varnames = colnames(Y)
    ),
    class = "lp_forecast"
  )
}

#' @export
print.lp_forecast <- function(x, digits = 4, ...) {
  cat("<lp_forecast> direct forecasts, horizons 1..", x$h, "\n", sep = "")
  print(round(x$mean, digits))
  invisible(x)
}

#' Compare local-projection and VAR impulse responses
#'
#' Plagborg-Moller & Wolf (2021): with the same lag length and the same recursive identification,
#' population LP and VAR responses coincide at all horizons; in finite samples they differ only by
#' estimation error. This helper estimates the LP \eqn{y_{t+h}=\alpha_h+\beta_hy_{j,t}+\gamma_h'w_t+e}
#' with controls the contemporaneous variables ordered before the shock \eqn{j} and \eqn{p} lags of
#' all variables, and the VAR(\eqn{p}) Cholesky response (ordering = column order) normalised to a
#' unit impact on the shock variable, \eqn{\Theta_h e_j/\Theta_0^{jj}}, and returns both with their
#' difference.
#'
#' **Economic intuition.** Choosing between LP and VAR is a bias-variance trade-off: the VAR
#' extrapolates the first \eqn{p} lags to all horizons (low variance, bias if the lag order is too
#' short for long-run dynamics), LP is unbiased for any \eqn{h} but has variance growing with
#' \eqn{h}. A divergence at long horizons that is large relative to `lp_se` signals VAR
#' misspecification (or LP noise); agreement supports either. *Pitfall:* the comparison is only
#' meaningful with identical lags and ordering, which this function enforces. VAR responses are
#' point estimates (no bands).
#'
#' @inheritParams lp
#' @param shock <[`tidy-select`][tidyselect::language]> the shock variable of `Y`.
#' @return a `lp_var_cmp` tibble with `horizon`, `response`, `lp`, `lp_se`, `var`, `diff = lp - var`,
#'   `z = diff / lp_se`.
#' @references Plagborg-Moller, M. and Wolf, C. (2021) Local projections and VARs estimate the same
#'   impulse responses. *Econometrica*.
#' @family local projections
#' @examples
#' compare_lp_var(yraw, shock = unemp, horizon = 12, lags = 2)
#' @export
compare_lp_var <- function(
  Y,
  shock,
  horizon = 20L,
  lags = 4L,
  se = c("nw", "hc"),
  nw_lag = NULL
) {
  se <- match.arg(se)
  Y <- prep_data(Y)$Y
  vn <- colnames(Y)
  sj <- select_vars(rlang::enquo(shock), vn, "shock")
  if (length(sj) != 1L) {
    rlang::abort("`shock` must select exactly one variable")
  }
  H <- as.integer(horizon)
  before <- if (sj > 1L) seq_len(sj - 1L)
  fit <- lp(Y, !!sj, H, lags, contemp = !!before, se = se, nw_lag = nw_lag)
  v <- irf(
    varm(Y, lags),
    id_chol(),
    horizon = H,
    ci = "none",
    normalise = list(shock = sj, variable = sj, size = 1)
  )
  out <- tibble::tibble(
    horizon = rep(0:H, times = length(vn)),
    response = rep(vn, each = H + 1L),
    lp = as.vector(fit$estimate),
    lp_se = as.vector(fit$se),
    var = as.vector(t(v$estimate[, sj, ]))
  )
  out$diff <- out$lp - out$var
  out$z <- ifelse(out$lp_se > 0, out$diff / out$lp_se, NA_real_)
  structure(
    out,
    shock = vn[sj],
    lags = lags,
    class = c("lp_var_cmp", class(out))
  )
}

#' @export
print.lp_var_cmp <- function(x, digits = 3, ...) {
  cat(
    "<lp_var_cmp> shock '",
    attr(x, "shock"),
    "', ",
    attr(x, "lags"),
    " lags; mean |LP - VAR| / se_LP by horizon block:\n",
    sep = ""
  )
  block <- cut(x$horizon, 3, labels = c("short", "medium", "long"))
  print(round(
    tapply(abs(x$z), list(block, x$response), mean, na.rm = TRUE),
    digits
  ))
  invisible(x)
}
