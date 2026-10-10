#' Prior for the TVP-VAR (Primiceri training-sample calibration)
#'
#' Hyperparameters of the training-sample prior of Primiceri (2005). Let \eqn{\hat B,\hat\alpha,\hat\sigma}
#' be OLS estimates on the first `tau` observations (after lag loss), with covariance matrices
#' \eqn{V(\hat B), V(\hat\alpha)}. Then
#' \deqn{\beta_0\sim N(\hat B,\, c\,V(\hat B)),\quad \alpha_0\sim N(\hat\alpha,\, c\,V(\hat\alpha)),\quad
#'   \log\sigma_{i0}\sim N(\log\hat\sigma_i, 10),}
#' \deqn{Q\sim IW(k_Q^2\,\nu_Q\,V(\hat B),\nu_Q),\quad S_i\sim IW(k_S^2\,(1+\dim S_i)\,V(\hat\alpha_i),1+\dim S_i),\quad
#'   W_i\sim IG(1, k_W^2).}
#' The `k` multipliers are the prior standard deviations of the *per-period* drift relative to the
#' OLS standard errors; they are the single most influential choice in a TVP-VAR.
#'
#' **Economic intuition and pitfalls.** `kQ` encodes how fast the transmission mechanism (lag
#' coefficients) is believed to move. With `kQ = 0.01` the prior says a coefficient drifts by about 1% of
#' its OLS standard error per quarter: over 170 quarters that is a standard deviation of roughly 13%
#' of the standard error times \eqn{\sqrt{170}}, so substantial but gradual change is allowed. Raising
#' `kQ` lets the data attribute any misfit (outliers, omitted variables, mis-specified lags) to
#' coefficient drift, i.e. *spurious* time variation. The prior is not neutral: with an
#' inverse-Wishart prior whose degrees of freedom are tied to a short training sample, the posterior of `Q`
#' is dominated by the prior in samples of typical macro length. Always report results for several
#' `kQ` values and compare with the constant-coefficient fit.
#'
#' @param tau training-sample length (effective observations, i.e. after the first `p`).
#' @param kQ,kS,kW prior scale multipliers for the innovations of the coefficients, of the
#'   contemporaneous (A) elements and of the log-volatilities.
#' @param vscale multiplier `c` on the OLS covariance in the prior for the initial states.
#' @param df_Q prior degrees of freedom of `Q` (default `tau`).
#' @param h_var prior variance of the initial log-variance (default 10).
#' @return a `tvp_prior` object for [tvpvar()].
#' @family TVP-VAR
#' @export
tvp_prior <- function(
  tau = 40L,
  kQ = 0.01,
  kS = 0.1,
  kW = 0.01,
  vscale = 4,
  df_Q = NULL,
  h_var = 10
) {
  structure(
    list(
      tau = as.integer(tau),
      kQ = kQ,
      kS = kS,
      kW = kW,
      vscale = vscale,
      df_Q = df_Q,
      h_var = h_var
    ),
    class = "tvp_prior"
  )
}

# OLS training-sample moments on the first `tau` observations
tvp_train <- function(Y, X, K, m, prior, q_diag) {
  tau <- prior$tau
  Yt <- Y[seq_len(tau), , drop = FALSE]
  Xt <- X[seq_len(tau), , drop = FALSE]
  XtXi <- solve(crossprod(Xt))
  B <- XtXi %*% crossprod(Xt, Yt)
  U <- Yt - Xt %*% B
  Sig <- crossprod(U) / (tau - m)
  VB <- kronecker(Sig, XtXi)
  a0 <- numeric(0)
  Pa <- vector("list", K)
  s2 <- numeric(K)
  s2[1] <- Sig[1, 1] * (tau - m) / tau * tau / (tau - m) # var of u_1
  for (i in seq_len(K)[-1]) {
    E <- U[, seq_len(i - 1), drop = FALSE]
    EtEi <- solve(crossprod(E))
    cf <- EtEi %*% crossprod(E, U[, i])
    v <- sum((U[, i] - E %*% cf)^2) / (tau - i + 1)
    a0 <- c(a0, -as.numeric(cf))
    Pa[[i]] <- v * EtEi
    s2[i] <- v
  }
  nuQ <- if (is.null(prior$df_Q)) tau else prior$df_Q
  Q0 <- prior$kQ^2 * VB
  if (q_diag) {
    Q0 <- diag(diag(Q0), nrow(Q0))
  }
  S0 <- Sbar <- vector("list", K)
  nuS <- numeric(K)
  for (i in seq_len(K)[-1]) {
    nuS[i] <- i # 1 + dim(S_i), dim = i - 1
    Sbar[[i]] <- prior$kS^2 * nuS[i] * Pa[[i]]
    S0[[i]] <- prior$kS^2 * Pa[[i]]
    Pa[[i]] <- prior$vscale * Pa[[i]]
  }
  list(
    b0 = as.numeric(B),
    P0 = prior$vscale * VB,
    a0 = a0,
    Pa0 = Pa,
    Qbar = prior$kQ^2 * nuQ * VB,
    nuQ = nuQ,
    Q0 = Q0,
    Sbar = Sbar,
    nuS = nuS,
    S0 = S0,
    hmu = log(s2)
  )
}

#' Time-varying parameter VAR with stochastic volatility
#'
#' Bayesian TVP-VAR of Cogley & Sargent (2001, 2005) and Primiceri (2005), with the Gibbs-sampler
#' ordering of Del Negro & Primiceri (2015):
#' \deqn{y_t = X_t'B_t + u_t,\quad u_t = A_t^{-1}\Lambda_t^{1/2}\varepsilon_t,\quad \varepsilon_t\sim N(0,I_K),}
#' \deqn{\beta_t=\mathrm{vec}(B_t)=\beta_{t-1}+\eta_t,\ \eta_t\sim N(0,Q);\quad
#'   \alpha_t=\alpha_{t-1}+\zeta_t,\ \zeta_t\sim N(0,S);\quad
#'   \log\lambda_{it}=\log\lambda_{i,t-1}+\nu_{it},\ \nu_{it}\sim N(0,W_i),}
#' where \eqn{X_t=(1,y_{t-1}',\dots,y_{t-p}')'}, \eqn{A_t} is unit lower triangular with free elements
#' \eqn{\alpha_t}, \eqn{\Sigma_t=A_t^{-1}\Lambda_tA_t^{-1\prime}} and \eqn{S} is block diagonal with one
#' block per row of \eqn{A_t}. Each Gibbs sweep draws, in order, \eqn{\beta^T|A^T,\Lambda^T,Q},
#' \eqn{\alpha^T|\beta^T,\Lambda^T,S}, the Kim-Shephard-Chib (Omori et al. 2007, 10-component) mixture
#' indicators and \eqn{\log\lambda^T|\beta^T,\alpha^T,W}, and finally \eqn{Q,S,W} from their (inverse-)
#' Wishart / gamma conditionals. State paths are drawn with the Durbin-Koopman (2002) simulation
#' smoother (cost linear in \eqn{T}, quadratic in the state dimension), all in C++ with R's RNG.
#'
#' Special cases: `tv_A = FALSE` gives the Cogley-Sargent (2005) model (time-varying coefficients and volatilities,
#' constant contemporaneous matrix); `sv = FALSE, tv_A = FALSE` is the Cogley-Sargent (2001) model with homoskedastic
#' innovations; `tv_beta = FALSE` freezes the coefficients (a Bayesian VAR with SV, the benchmark for judging
#' whether time variation is needed). The recursive ordering of the variables matters whenever `tv_A` or
#' `sv` is on, because \eqn{A_t} imposes a Cholesky structure on \eqn{\Sigma_t}.
#'
#' **Economic intuition.** Time variation in \eqn{B_t} is a changing *transmission mechanism* (for example a
#' monetary policy rule that responds more strongly to inflation after Volcker); time variation in
#' \eqn{\Lambda_t} is a change in the size of shocks (the Great Moderation, the 2008 and 2020 episodes); time variation
#' in \eqn{A_t} is a change in the contemporaneous propagation of shocks. These are conceptually distinct, and a
#' model that omits stochastic volatility will tend to push heteroskedasticity into the coefficients. **Pitfalls:**
#' (i) the prior on `Q` (see [tvp_prior()]) is informative and drives how much coefficient variation is found;
#' (ii) with many lags the state dimension (\eqn{K(Kp+1)}) is large relative to the sample, so
#' keep `p` small; (iii) the first `tau` observations are used up as a training sample and are not part
#' of the estimation sample; (iv) draws of \eqn{B_t} may be explosive at some dates: use
#' `reject_explosive = TRUE` for the (slow) Cogley-Sargent stability truncation.
#'
#' @inheritParams varm
#' @param vars <[`tidy-select`][tidyselect::language]> endogenous variables
#'   (default: all numeric columns). A `Date` column of a data frame supplies the
#'   dates of the paths.
#' @param ndraw number of retained draws.
#' @param nburn burn-in iterations.
#' @param thin keep every `thin`-th sweep.
#' @param sv stochastic volatility (`TRUE`) or constant diagonal \eqn{\Lambda} (`FALSE`).
#' @param tv_beta,tv_A time-varying coefficients / contemporaneous matrix \eqn{A_t}.
#' @param q_diag diagonal \eqn{Q} (independent inverse-gamma priors) instead of a full inverse-Wishart.
#' @param prior a [tvp_prior()].
#' @param nsave_paths number of (evenly spaced) retained draws for which the complete state paths are stored
#'   (memory is `nsave_paths x T x K(Kp+1)` doubles); impulse responses, volatilities and forecasts are computed
#'   from these. Posterior means of the paths use all draws.
#' @param reject_explosive redraw the coefficient path (up to `max_tries` times, otherwise keep the previous
#'   path) when any date has a companion root outside the unit circle.
#' @param max_tries see `reject_explosive`.
#' @return a `tvpvar` object. Methods: [print()], [summary()], [coef()], [irf()]
#'   (`date=`), [fevd()], [forecast()], [tidy()][tidy.tvpvar()],
#'   [autoplot()][autoplot.tvpvar()], plus [irf_path()] and [tvp_volatility()].
#' @references Cogley & Sargent (2005), Primiceri (2005), Del Negro & Primiceri (2015), Durbin & Koopman (2002).
#' @family TVP-VAR
#' @examples
#' \donttest{
#' fit <- tvpvar(yraw, p = 2, ndraw = 500, nburn = 500)
#' fit
#' autoplot(fit, "vol")
#' autoplot(fit, "irf", shock = rate, response = infl, horizon = c(0, 4, 8))
#' plot(irf(fit, date = 50, horizon = 12))
#' }
#' @export
tvpvar <- function(
  Y,
  p = 2L,
  ndraw = 5000L,
  nburn = 2000L,
  thin = 1L,
  const = TRUE,
  sv = TRUE,
  tv_beta = TRUE,
  tv_A = TRUE,
  q_diag = FALSE,
  prior = tvp_prior(),
  nsave_paths = min(ndraw, 1000L),
  reject_explosive = FALSE,
  max_tries = 10L,
  vars = NULL
) {
  cl <- match.call()
  if (!inherits(prior, "tvp_prior")) {
    rlang::abort("`prior` must be a tvp_prior()")
  }
  ts_time <- if (stats::is.ts(Y)) as.numeric(stats::time(Y))
  input <- prep_data(Y, rlang::enquo(vars))
  Y <- input$Y
  p <- as.integer(p)
  K <- ncol(Y)
  Tn <- nrow(Y)
  det <- make_det(Tn, const, FALSE, NULL)
  ols <- cpp_var_ols(Y, p, det)
  Ye <- ols$Y
  X <- ols$X
  m <- ncol(X)
  tau <- prior$tau
  if (tau <= m + K) {
    rlang::abort(paste0(
      "training sample `tau` must exceed the number of regressors per equation + K (",
      m + K,
      ")"
    ))
  }
  if (nrow(Ye) - tau < 20) {
    rlang::abort("estimation sample too short after the training sample")
  }
  pr <- tvp_train(Ye, X, K, m, prior, q_diag)
  if (pr$nuQ <= m * K + 1) {
    rlang::warn(paste0(
      "prior degrees of freedom for Q (",
      pr$nuQ,
      ") do not exceed the number of coefficients + 1 (",
      m * K + 1,
      "): the prior mean of Q is undefined; raise `tau` or `df_Q`"
    ))
  }
  est <- (tau + 1L):nrow(Ye)
  all_dates <- input$dates
  if (is.null(all_dates)) {
    all_dates <- if (is.null(ts_time)) seq_len(Tn) else ts_time
  }
  dates <- all_dates[p + est]
  vn <- colnames(Y)
  rn <- coef_names(vn, p, colnames(det))
  opts <- c(sv, tv_beta, tv_A, q_diag, reject_explosive, max_tries)
  out <- cpp_tvp(
    Ye[est, , drop = FALSE],
    X[est, , drop = FALSE],
    p,
    pr$b0,
    pr$P0,
    pr$a0,
    pr$Pa0,
    pr$Qbar,
    pr$nuQ,
    pr$Q0,
    pr$Sbar,
    pr$nuS,
    pr$S0,
    pr$hmu,
    prior$h_var,
    c(1, prior$kW^2),
    c(3, 2),
    as.numeric(opts),
    as.integer(ndraw),
    as.integer(nburn),
    as.integer(thin),
    as.integer(nsave_paths)
  )
  Te <- length(est)
  nsp <- dim(out$beta)[3]
  out$beta <- array(
    out$beta,
    c(m, K, Te, nsp),
    dimnames = list(rn, vn, NULL, NULL)
  )
  out$beta_mean <- array(
    out$beta_mean,
    c(m, K, Te),
    dimnames = list(rn, vn, NULL)
  )
  info <- list(
    Y = Y,
    p = p,
    K = K,
    m = m,
    T = Te,
    varnames = vn,
    coefnames = rn,
    dates = dates,
    const = const,
    sv = sv,
    tv_beta = tv_beta,
    tv_A = tv_A,
    q_diag = q_diag,
    prior = prior,
    prior_values = pr,
    ndraw = ndraw,
    tau = tau,
    call = cl
  )
  structure(c(out, info), class = "tvpvar")
}

#' @export
print.tvpvar <- function(x, ...) {
  cat(
    "<tvpvar> TVP-VAR(",
    x$p,
    "), K = ",
    x$K,
    ", estimation sample T = ",
    x$T,
    " (",
    format(x$dates[1]),
    " - ",
    format(x$dates[x$T]),
    "), training tau = ",
    x$tau,
    "\n",
    sep = ""
  )
  cat(
    "  coefficients: ",
    if (x$tv_beta) "time-varying" else "constant",
    "; A_t: ",
    if (x$tv_A) "time-varying" else "constant",
    "; volatility: ",
    if (x$sv) "stochastic" else "constant",
    "\n",
    sep = ""
  )
  cat(
    "  ",
    x$ndraw,
    " retained draws, ",
    dim(x$h)[3],
    " stored paths\n",
    sep = ""
  )
  invisible(x)
}

#' @export
summary.tvpvar <- function(object, ...) {
  nsp <- dim(object$h)[3]
  sds <- cpp_tvp_sd(object$alpha, object$h, object$K)[,,
    seq_len(nsp),
    drop = FALSE
  ]
  med_range <- function(a) range(apply(a, 1, stats::median))
  structure(
    list(
      object = object,
      q_sd = if (object$tv_beta) {
        stats::median(rowMeans(object$q_sd))
      } else {
        NA_real_
      },
      w_sd = stats::setNames(colMeans(object$w_sd), object$varnames),
      s_sd = colMeans(object$s_sd),
      sd_range = t(apply(sds, 1, med_range))
    ),
    class = "summary.tvpvar"
  )
}

#' @export
print.summary.tvpvar <- function(x, ...) {
  print(x$object)
  cat("Posterior mean innovation s.d. of the log-volatilities:\n")
  print(round(x$w_sd, 4))
  if (!is.na(x$q_sd)) {
    cat(
      "Median (over coefficients) posterior-mean s.d. of the coefficient drift sqrt(Q_jj): ",
      format(x$q_sd, digits = 3),
      "\n",
      sep = ""
    )
  }
  dimnames(x$sd_range) <- list(x$object$varnames, c("min", "max"))
  cat("Range over time of the median reduced-form residual s.d.:\n")
  print(round(x$sd_range, 3))
  invisible(x)
}

# Positions of `date` (nearest match, or "first"/"last") in `object$dates`
tvp_date_index <- function(object, date) {
  if (is.null(date)) {
    return(object$T)
  }
  if (is.character(date) && all(date %in% c("first", "last"))) {
    return(ifelse(date == "last", object$T, 1L))
  }
  if (inherits(object$dates, "Date")) {
    date <- as.Date(date)
  }
  if (any(date < min(object$dates) | date > max(object$dates))) {
    rlang::warn(
      "`date` lies outside the estimation sample; using the nearest date"
    )
  }
  vapply(
    date,
    function(d) which.min(abs(object$dates - d)),
    1L,
    USE.NAMES = FALSE
  )
}

#' @export
coef.tvpvar <- function(object, date = NULL, stat = c("mean", "median"), ...) {
  stat <- match.arg(stat)
  a <- if (stat == "mean") {
    object$beta_mean
  } else {
    apply(object$beta, 1:3, stats::median)
  }
  if (is.null(date)) {
    return(a)
  }
  idx <- tvp_date_index(object, date)
  if (length(idx) == 1L) a[,, idx] else a[,, idx, drop = FALSE]
}

#' Time-varying impulse responses, variance decompositions and forecasts of a TVP-VAR
#'
#' Impulse responses are computed *at a date* by freezing the parameters at their date-`t` values
#' (the standard "local" IRF of Cogley-Sargent / Primiceri):
#' \deqn{\Theta_h(t)=\Phi_h(B_t)\,A_t^{-1}\Lambda_t^{1/2},\qquad \Phi_h=\sum_{l=1}^{\min(h,p)}A_l(t)\Phi_{h-l},}
#' i.e. the response of \eqn{y_{t+h}} to a one-standard-deviation recursive shock if the economy were to stay at the
#' date-\eqn{t} parameters. Note that the shock *size* \eqn{\Lambda_t^{1/2}} varies with the date, so changes in the
#' IRF mix changes in propagation and in volatility; use `normalise = list(shock =, variable =, size =)` to
#' hold the impact response fixed (e.g. a 25bp rate shock) and isolate changes in propagation. The
#' FEVD uses the same frozen-parameter responses. The identification is recursive in the order of the data.
#'
#' Forecasts simulate, for every stored posterior draw, future paths of \eqn{\beta,\alpha,\log\lambda} as random
#' walks with the draw's own \eqn{Q,S,W} (`drift = TRUE`), so predictive densities include the possibility that
#' parameters keep moving; `drift = FALSE` freezes them at the end-of-sample values.
#'
#' **Economic intuition.** Comparing `irf(fit, date = c(1975, 2005))` asks whether, say, the output response to a
#' monetary tightening was larger in the 1970s than in the Great Moderation. Bands reflect posterior uncertainty
#' about the date-specific parameters only; they do not reflect uncertainty about how long the parameters stay put.
#'
#' @param object a [tvpvar()] fit.
#' @param date one or several dates (matched to the nearest element of `object$dates`: the `ts` time stamp
#'   if `Y` was a `ts`, else the row index of `Y`), or `"last"`/`"first"`. Default: the last date.
#' @inheritParams dynamics
#' @param ... unused.
#' @return for `irf`/`fevd`: a `varm_irf`/`varm_fevd` (single date, plottable) or a named list of them of class
#'   `tvpvar_irf` (several dates). `forecast` returns a `varm_forecast`.
#' @family TVP-VAR
#' @name tvp_dynamics
NULL

#' @rdname tvp_dynamics
#' @export
irf.tvpvar <- function(
  object,
  date = NULL,
  horizon = 20L,
  level = c(0.68, 0.90),
  normalise = NULL,
  cumulative = NULL,
  ...,
  .fevd = FALSE
) {
  if ("id" %in% names(list(...))) {
    rlang::abort(
      "`tvpvar` IRFs are recursive (Cholesky) only; `id` is not supported"
    )
  }
  K <- object$K
  H <- as.integer(horizon)
  idx <- tvp_date_index(object, date)
  N <- dim(object$h)[3]
  a <- cpp_tvp_irf(
    array(object$beta, c(object$m * K, object$T, N)),
    object$alpha,
    object$h,
    object$m,
    K,
    object$p,
    H,
    idx - 1L
  )
  a <- array(a, c(K, K, H + 1L, N, length(idx)))
  one <- function(k) {
    dyn_object(
      a[,,,, k],
      K,
      object$p,
      H,
      id_chol(),
      object$varnames,
      level,
      normalise,
      cumulative,
      fevd = .fevd,
      extra = list(date = object$dates[idx[k]], method = "tvpvar")
    )
  }
  res <- lapply(seq_along(idx), one)
  if (length(idx) == 1L) {
    return(res[[1]])
  }
  names(res) <- format(object$dates[idx])
  structure(res, class = "tvpvar_irf")
}

#' @rdname tvp_dynamics
#' @export
fevd.tvpvar <- function(
  object,
  date = NULL,
  horizon = 20L,
  level = c(0.68, 0.90),
  ...
) {
  irf.tvpvar(object, date, horizon, level, ..., .fevd = TRUE)
}

#' @export
print.tvpvar_irf <- function(x, ...) {
  cat(
    "<tvpvar_irf> responses at ",
    length(x),
    " dates: ",
    paste(names(x), collapse = ", "),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' @rdname tvp_dynamics
#' @param h forecast horizon (steps).
#' @param drift simulate future parameter drift (`TRUE`) or freeze parameters at the last date.
#' @export
forecast.tvpvar <- function(
  object,
  h = 8L,
  level = c(0.68, 0.90),
  drift = TRUE,
  ...
) {
  K <- object$K
  p <- object$p
  Tn <- nrow(object$Y)
  h <- as.integer(h)
  Tl <- object$T
  f <- cpp_tvp_forecast(
    matrix(object$beta[,, Tl, ], object$m * K),
    matrix(object$alpha[, Tl, ], nrow(object$alpha)),
    matrix(object$h[, Tl, ], K),
    object$Q,
    object$S,
    object$W,
    object$m,
    K,
    p,
    object$Y[(Tn - p + 1):Tn, , drop = FALSE],
    matrix(1, h, as.integer(object$const)),
    h,
    drift,
    object$tv_beta,
    object$tv_A,
    object$sv
  )
  probs <- level_probs(level)
  bands <- aperm(apply(f, c(2, 3), stats::quantile, probs = probs), c(2, 3, 1))
  dimnames(bands) <- list(NULL, object$varnames, format(probs))
  fmean <- apply(f, c(2, 3), mean)
  colnames(fmean) <- object$varnames
  structure(
    list(
      mean = fmean,
      bands = bands,
      draws = f,
      h = h,
      varnames = object$varnames,
      level = level
    ),
    class = "varm_forecast"
  )
}

# Posterior median and equal-tailed band of draws along the last dimension of `a`
# (array); returns array (3, dim(a)[-last]) of lower, median, upper
tvp_quantiles <- function(a, level) {
  k <- length(dim(a))
  pr <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  apply(a, seq_len(k - 1L), stats::quantile, probs = pr)
}

#' Response of one variable to one shock at a fixed horizon, across all dates
#'
#' Evolution over time of \eqn{\Theta_{h}^{ij}(t)}, the response of variable \eqn{i} to recursive shock \eqn{j} at
#' horizon \eqn{h} computed with frozen date-\eqn{t} parameters (see [irf.tvpvar()]). This is the usual "time-varying
#' IRF" plot: for example the impact (or 8-quarter) response of inflation to a policy shock through time. With
#' `normalise_to` the shock is rescaled each date so the impact response of `variable` equals `size`
#' (e.g. 1 percentage point on the policy rate), removing changes in shock size.
#'
#' @param object a [tvpvar()] fit.
#' @param shock,response <[`tidy-select`][tidyselect::language]> one variable
#'   each.
#' @param horizon horizon of the response (0 = impact).
#' @param normalise_to optional <[`tidy-select`][tidyselect::language]> variable
#'   whose impact response is normalised to `size`.
#' @param size normalised impact (default 1).
#' @param level credible level of the band.
#' @return a `tvpvar_path` tibble with columns `date`, `estimate` (posterior
#'   median), `conf.low`, `conf.high`; [autoplot()][autoplot.tvpvar()] plots it.
#' @family TVP-VAR
#' @examples
#' \donttest{
#' fit <- tvpvar(yraw, p = 2, ndraw = 200, nburn = 200)
#' irf_path(fit, shock = rate, response = infl, horizon = 4)
#' }
#' @export
irf_path <- function(
  object,
  shock,
  response,
  horizon = 0L,
  normalise_to = NULL,
  size = 1,
  level = 0.68
) {
  if (!inherits(object, "tvpvar")) {
    rlang::abort("`object` must be a tvpvar fit")
  }
  vn <- object$varnames
  one <- function(quo, arg) {
    j <- select_vars(quo, vn, arg)
    if (length(j) != 1L) {
      rlang::abort(paste0("`", arg, "` must select exactly one variable"))
    }
    j
  }
  j <- one(rlang::enquo(shock), "shock")
  i <- one(rlang::enquo(response), "response")
  nv <- select_or(rlang::enquo(normalise_to), vn, "normalise_to", NULL)
  r <- cpp_tvp_irf_path(
    array(object$beta, c(object$m * object$K, object$T, dim(object$h)[3])),
    object$alpha,
    object$h,
    object$m,
    object$K,
    object$p,
    as.integer(horizon),
    i - 1L,
    j - 1L,
    if (is.null(nv)) -1L else nv - 1L,
    size
  )
  q <- tvp_quantiles(r, level)
  out <- tibble::tibble(
    date = object$dates,
    estimate = q[2, ],
    conf.low = q[1, ],
    conf.high = q[3, ]
  )
  structure(
    out,
    label = paste0(vn[i], " <- ", vn[j], ", h = ", horizon),
    class = c("tvpvar_path", class(out))
  )
}

#' Volatility paths of a TVP-VAR
#'
#' Posterior bands for the time-varying standard deviations: `type = "reduced"` gives
#' \eqn{\sqrt{\Sigma_{t,ii}}} (reduced-form residual s.d. of each variable, what a forecaster faces);
#' `"structural"` gives \eqn{\exp(h_{it}/2)}, the s.d. of the orthogonalised (recursive) shocks.
#' Comparing the two separates changes in shock size from changes in contemporaneous propagation \eqn{A_t}.
#'
#' @param object a [tvpvar()] fit.
#' @param type `"reduced"` or `"structural"`.
#' @param level credible level.
#' @return a tibble with `date`, `variable`, `estimate` (posterior median),
#'   `conf.low`, `conf.high`.
#' @family TVP-VAR
#' @export
tvp_volatility <- function(
  object,
  type = c("reduced", "structural"),
  level = 0.68
) {
  type <- match.arg(type)
  K <- object$K
  a <- if (type == "reduced") {
    cpp_tvp_sd(object$alpha, object$h, K)
  } else {
    exp(0.5 * object$h)
  }
  q <- tvp_quantiles(a, level) # 3 x K x T
  tibble::tibble(
    date = rep(object$dates, K),
    variable = factor(rep(object$varnames, each = object$T), object$varnames),
    estimate = as.vector(t(q[2, , ])),
    conf.low = as.vector(t(q[1, , ])),
    conf.high = as.vector(t(q[3, , ]))
  )
}
