#' State-dependent local projections (Auerbach-Gorodnichenko)
#'
#' Local projections whose coefficients depend on the state of the economy through a transition
#' function of a switching variable \eqn{z}. For each horizon \eqn{h},
#' \deqn{y_{t+h} = F(z_{t-1})\,[\alpha^L_h + \beta^L_h s_t + \gamma^{L\prime}_h w_t] +
#'   (1-F(z_{t-1}))\,[\alpha^H_h + \beta^H_h s_t + \gamma^{H\prime}_h w_t] + e_{t+h},}
#' with the smooth logistic transition
#' \eqn{F(z) = 1/(1+\exp\{\gamma(\tilde z - c)\})} on the standardised variable
#' \eqn{\tilde z=(z-\bar z)/\mathrm{sd}(z)} (Auerbach & Gorodnichenko 2012, 2013;
#' \eqn{F\to1} in the "low" state, e.g. recession when \eqn{z} is growth) or the binary version
#' \eqn{F = 1\{\tilde z_{t-1} < c\}}. The state enters lagged by one period to avoid contemporaneous
#' feedback of the shock on the state. All \eqn{2\times} regressors are estimated in one
#' regression per horizon (C++ engine), with Newey-West errors, and the cross-horizon covariance
#' is used to test \eqn{H_0:\beta^L_h=\beta^H_h} pointwise (HAC t-test) and jointly across
#' horizons (Wald).
#'
#' **Timing conventions.** `timing = "contemporaneous"` (default) is the shock-regression form above
#' (observed or external shock \eqn{s_t}, as in [lp()]). `timing = "lagged"` replicates the convention of
#' `lpirfs::lp_nl`: the regressors are lags of \eqn{y} interacted with the state, with no
#' contemporaneous terms; the horizon-\eqn{h} response is the coefficient on \eqn{y_{t-1}} in the
#' equation for \eqn{y_{t+h-1}} multiplied by column `shock` of the Cholesky impact matrix of a linear
#' VAR (standard-deviation shock when `shock_type = 0`, unit own-impact when `1`), and the impact response
#' (h = 0) is that column of the impact matrix in both states. Differences to `lpirfs`: bands here are
#' normal quantiles at `level` (not a user-supplied multiplier); `lpirfs` uses Newey-West lag
#' \eqn{h} for its horizon \eqn{h} (regression of \eqn{y_{t+h-1}}), which equals the default
#' \eqn{L=(h-1)+1} here; `lpirfs` has no state-specific constants (`state_const = FALSE`).
#'
#' **Economic intuition.** Fiscal and monetary multipliers are often larger in recessions; this
#' estimates state-specific responses without a nonlinear VAR, so the state can change over the
#' horizon without being imposed. The cost is variance: each state uses only the observations
#' weighted by \eqn{F} or \eqn{1-F}, so effective samples are small and bands wide, and with persistent
#' states the difference test has low power at long horizons. *Pitfalls:* the calibrated
#' `gamma` and cutoff are not estimated (AG use \eqn{\gamma=1.5} and a cutoff so that the low state has
#' about 20% weight, see `low_share`); a state that is endogenous to the shock beyond one lag
#' contaminates the comparison; with the binary version one state may have very few observations.
#'
#' @inheritParams lp
#' @param state switching variable \eqn{z}: a variable of `Y` (tidyselect, one
#'   variable) or a numeric vector of length T.
#' @param transition `"logistic"` (smooth) or `"binary"`.
#' @param gamma logistic smoothness (larger is closer to a regime switch).
#' @param cutoff cutoff \eqn{c} on the standardised scale; ignored if
#'   `low_share` is given.
#' @param low_share optional share in (0, 1) of observations to place below the
#'   cutoff (sets \eqn{c} to that quantile of \eqn{\tilde z}).
#' @param hp_lambda if not `NULL`, use the Hodrick-Prescott cycle of `state`
#'   (smoothing parameter, e.g. 1600 for quarterly data) as switching variable.
#' @param standardise standardise \eqn{z} (after filtering) to mean 0, sd 1.
#' @param timing `"contemporaneous"` or `"lagged"` (the `lpirfs` convention),
#'   see Details.
#' @param lags_lin lag length of the linear VAR used for the impact matrix
#'   (`timing = "lagged"`).
#' @param shock_type `"lagged"` only: `0` standard-deviation, `1` unit shock.
#' @param state_const interact the constant with the state (state-specific
#'   intercepts); `lpirfs` uses a common constant, so set `FALSE` when comparing.
#' @param trend `"lagged"` only: include a linear time trend (`lpirfs` trend = 1).
#' @return an `lp_state` object: arrays `$estimate`, `$se` (`(H+1) x K x 3`, last
#'   dimension `high`, `low`, `diff = low - high`), `$bands`, `$tstat`/`$pvalue`
#'   for the difference, `$wald` (joint test of equality over horizons, per
#'   response), the transition `$F` and `$cutoff`, and `$raw` (coefficient array
#'   and cross-horizon covariance from the engine) with `$contrasts`.
#'   [tidy()][tidy.lp_state()] returns one row per horizon, response, state and
#'   band level.
#' @references Auerbach, A. and Gorodnichenko, Y. (2012) Measuring the output
#'   responses to fiscal policy. *AEJ: Economic Policy*; (2013) Fiscal
#'   multipliers in recession and expansion.
#' @family local projections
#' @examples
#' fit <- lp_state(
#'   yraw,
#'   shock = rate,
#'   state = infl,
#'   horizon = 8,
#'   lags = 2,
#'   contemp = c(infl, unemp),
#'   gamma = 3,
#'   low_share = 0.3
#' )
#' fit$estimate[, 3, ]
#' @export
lp_state <- function(
  Y,
  shock,
  state,
  horizon = 12L,
  lags = 4L,
  contemp = NULL,
  controls = NULL,
  responses = NULL,
  cumulative = NULL,
  transition = c("logistic", "binary"),
  gamma = 1.5,
  cutoff = 0,
  low_share = NULL,
  hp_lambda = NULL,
  standardise = TRUE,
  timing = c("contemporaneous", "lagged"),
  lags_lin = lags,
  shock_type = 0L,
  state_const = TRUE,
  trend = FALSE,
  se = c("nw", "hc"),
  nw_lag = NULL,
  level = c(0.68, 0.90),
  dof_adjust = FALSE
) {
  transition <- match.arg(transition)
  timing <- match.arg(timing)
  se <- match.arg(se)
  Y <- prep_data(Y)$Y
  n <- nrow(Y)
  vn <- colnames(Y)
  H <- as.integer(horizon)
  lags <- as.integer(lags)
  if (lags < 1L) {
    rlang::abort("`lags` must be >= 1")
  }
  s <- lp_series(rlang::enquo(shock), Y, "shock")
  z <- lp_series(rlang::enquo(state), Y, "state")$x
  rs <- select_or(rlang::enquo(responses), vn, "response", seq_along(vn))
  cum <- as.integer(
    rs %in% select_or(rlang::enquo(cumulative), vn, "cumulative", integer())
  )
  contemp <- select_or(rlang::enquo(contemp), vn, "contemp", NULL)

  tr <- state_transition(
    z,
    transition,
    gamma,
    cutoff,
    low_share,
    hp_lambda,
    standardise
  )
  Fl <- c(NA_real_, tr$F[-n])
  if (timing == "lagged" && (!is.null(contemp) || !is.null(controls))) {
    rlang::abort(
      "`contemp`/`controls` are not supported with timing = 'lagged'"
    )
  }
  des <- if (timing == "lagged") {
    lp_state_lagged(Y, s, Fl, lags, lags_lin, shock_type, trend, rs, H)
  } else {
    lp_state_contemp(Y, s, Fl, lags, contemp, controls, state_const, H)
  }

  seT <- if (se == "nw") 0L else 1L
  nwf <- if (is.null(nw_lag)) -1L else as.integer(nw_lag)
  Yr <- Y[, rs, drop = FALSE]
  Z0 <- matrix(0, n, 0L)
  r <- cpp_lp_cov(Yr, des$S, des$W, Z0, des$Hm, lags, cum, seT, nwf, dof_adjust)
  dn <- list(0:H, vn[rs], c("high", "low", "diff"))
  paths <- state_paths(r, des, length(rs), H, dn)
  est <- paths$est
  sev <- paths$se

  # joint Wald test over horizons: common bandwidth for a well-behaved matrix
  rw <- r
  if (is.null(nw_lag)) {
    rw <- cpp_lp_cov(
      Yr,
      des$S,
      des$W,
      Z0,
      des$Hm,
      lags,
      cum,
      seT,
      as.integer(des$Hm + 1L),
      dof_adjust
    )
  }
  wald <- state_wald(rw, est, des, vn[rs])
  tstat <- matrix(est[,, 3] / sev[,, 3], H + 1L, length(rs), dimnames = dn[1:2])
  probs <- level_probs(level)
  structure(
    list(
      estimate = est,
      se = sev,
      bands = band_array(est, sev, probs),
      probs = probs,
      level = level,
      tstat = tstat,
      pvalue = 2 * stats::pnorm(-abs(tstat)),
      wald = wald,
      F = Fl,
      cutoff = tr$cutoff,
      gamma = gamma,
      transition = transition,
      share_low = mean(Fl > 0.5, na.rm = TRUE),
      timing = timing,
      raw = r,
      contrasts = des$Cm,
      nobs = r$nobs,
      H = H,
      shock = des$sname,
      responses = vn[rs],
      lags = lags
    ),
    class = "lp_state"
  )
}

# Transition function F(z) of the (filtered, standardised) switching variable
state_transition <- function(
  z,
  transition,
  gamma,
  cutoff,
  low_share,
  hp_lambda,
  standardise
) {
  if (!is.null(hp_lambda)) {
    z <- z - hp_trend(z, hp_lambda)
  }
  if (standardise) {
    z <- (z - mean(z)) / stats::sd(z)
  }
  if (!is.null(low_share)) {
    cutoff <- unname(stats::quantile(z, low_share))
  }
  Fz <- if (transition == "logistic") {
    1 / (1 + exp(gamma * (z - cutoff)))
  } else {
    as.numeric(z < cutoff)
  }
  list(F = Fz, cutoff = cutoff)
}

hp_trend <- function(z, lambda) {
  n <- length(z)
  D <- diff(diag(n), differences = 2)
  as.numeric(solve(diag(n) + lambda * crossprod(D), z))
}

# Shock-regression form: shock and controls interacted with F(z_{t-1})
lp_state_contemp <- function(
  Y,
  s,
  Fl,
  lags,
  contemp,
  controls,
  state_const,
  H
) {
  d <- lp_design(Y, s, lags, contemp, controls)
  Wc <- d$W[, -1, drop = FALSE]
  W <- if (state_const) {
    cbind(d$W[, 1] * (1 - Fl), d$W[, 1] * Fl, Wc * (1 - Fl), Wc * Fl)
  } else {
    cbind(1, Wc * (1 - Fl), Wc * Fl)
  }
  list(
    S = cbind(d$s * (1 - Fl), d$s * Fl),
    W = W,
    Hm = H,
    off = 0L,
    imp = NULL,
    Cm = cbind(high = c(1, 0), low = c(0, 1), diff = c(-1, 1)),
    sname = d$sname
  )
}

# lpirfs::lp_nl convention: lags of Y interacted with the state, responses scaled
# by the Cholesky impact column of a linear VAR
lp_state_lagged <- function(
  Y,
  s,
  Fl,
  lags,
  lags_lin,
  shock_type,
  trend,
  rs,
  H
) {
  if (is.null(s$pos)) {
    rlang::abort("timing = 'lagged' needs a `shock` variable of `Y`")
  }
  if (H < 1L) {
    rlang::abort("timing = 'lagged' needs horizon >= 1")
  }
  n <- nrow(Y)
  K0 <- ncol(Y)
  L1 <- lag_matrix(Y, lags)
  S1 <- L1[, seq_len(K0), drop = FALSE]
  Wr <- if (lags > 1L) L1[, -seq_len(K0), drop = FALSE] else matrix(0, n, 0L)
  tr <- if (trend) seq_len(n)
  pl <- as.integer(lags_lin)
  Xl <- cbind(1, tr, lag_matrix(Y, pl))[(pl + 1):n, , drop = FALSE]
  U <- stats::lm.fit(Xl, Y[(pl + 1):n, , drop = FALSE])$residuals
  A <- t(chol(stats::cov(U)))
  j <- s$pos
  scale_j <- if (shock_type == 0) sqrt(diag(stats::cov(U)))[j] else 1
  impf <- A[, j] / A[j, j] * scale_j
  # contrasts act on the 2K lag-1 coefficients of equation k: the full impact vector
  Cm <- cbind(high = c(impf, rep(0, K0)), low = c(rep(0, K0), impf))
  list(
    S = cbind(S1 * (1 - Fl), S1 * Fl),
    W = cbind(1, tr, Wr * (1 - Fl), Wr * Fl),
    Hm = H - 1L,
    off = 1L,
    imp = impf[rs],
    Cm = cbind(Cm, diff = Cm[, 2] - Cm[, 1]),
    sname = colnames(Y)[j]
  )
}

# State-specific responses and standard errors: contrasts of the engine output
state_paths <- function(r, des, K, H, dn) {
  est <- se <- array(NA_real_, c(H + 1L, K, 3), dimnames = dn)
  nh <- des$Hm + 1L
  ns <- ncol(des$S)
  Cm <- des$Cm
  for (k in seq_len(K)) {
    for (hh in seq_len(nh)) {
      ii <- (hh - 1L) * ns + seq_len(ns)
      V <- r$vcov[ii, ii, k]
      est[des$off + hh, k, ] <- drop(crossprod(Cm, r$beta[hh, k, ]))
      se[des$off + hh, k, ] <- sqrt(pmax(diag(crossprod(Cm, V %*% Cm)), 0))
    }
  }
  if (!is.null(des$imp)) {
    for (k in seq_len(K)) {
      est[1, k, ] <- c(des$imp[k], des$imp[k], 0)
      se[1, k, ] <- 0
    }
  }
  list(est = est, se = se)
}

state_wald <- function(rw, est, des, rn) {
  nh <- des$Hm + 1L
  hidx <- des$off + seq_len(nh)
  Cd <- kronecker(diag(nh), des$Cm[, 3, drop = FALSE])
  wald <- vapply(
    seq_along(rn),
    function(k) wald_stat(est[hidx, k, 3], crossprod(Cd, rw$vcov[,, k] %*% Cd)),
    numeric(3)
  )
  wald <- t(wald)
  rownames(wald) <- rn
  wald
}

#' @export
print.lp_state <- function(x, digits = 3, ...) {
  cat(
    "<lp_state> ",
    x$timing,
    " timing, ",
    x$transition,
    " transition, shock '",
    x$shock,
    "', horizons 0..",
    x$H,
    "; low-state share ",
    round(x$share_low, 2),
    "\n",
    sep = ""
  )
  cat("Joint test of equal responses across states (Wald):\n")
  print(round(x$wald, digits))
  invisible(x)
}
