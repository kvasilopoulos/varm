#' Minnesota / sum-of-coefficients / single-unit-root prior
#'
#' Prior of Giannone, Lenza & Primiceri (2015) for a conjugate Normal-inverse-Wishart BVAR:
#' \deqn{\beta|\Sigma \sim N(b,\Sigma\otimes\Omega),\quad \Sigma\sim IW(\Psi, M+2).}
#' The prior variance of the coefficient on lag \eqn{l} of variable \eqn{j} in equation \eqn{i}
#' is \eqn{\lambda^2/(l^\alpha\psi_j)\,\psi_i}-scaled through \eqn{\Sigma}; the constant is
#' (nearly) flat. Dummy-observation priors shrink towards unit roots (`soc`, Doan-Litterman-Sims)
#' and towards cointegration/no-trend (`sur`, Sims).
#'
#' Each hyperparameter is a `c(mode, sd, min, max)` vector (gamma hyperprior with that mode and
#' standard deviation; hard bounds). `hyper` names the hyperparameters estimated hierarchically;
#' all others stay fixed at their `mode`.
#'
#' @param lambda overall tightness (`mode` 0.2).
#' @param alpha lag-decay exponent (`mode` 2; fixed unless in `hyper`).
#' @param soc sum-of-coefficients tightness \eqn{\mu} (smaller = tighter), or `NULL` to drop.
#' @param sur single-unit-root / dummy-initial-observation tightness \eqn{\theta}, or `NULL` to drop.
#' @param psi `NULL` (use AR(p) residual variances) or a numeric vector of length K
#'   (the scale of the inverse-Wishart prior, variance units).
#' @param delta prior mean of the own first lag: `1` (random walk, levels) or `0`
#'   (stationary, growth rates); scalar or length K.
#' @param hyper character vector among `"lambda"`, `"alpha"`, `"soc"`, `"sur"`, `"psi"`; or `"full"`.
#'   Default `c("lambda", "soc", "sur")`; `character(0)` fixes everything.
#' @param var prior variance of deterministic-term coefficients.
#' @return a `prior_minnesota` object for [bvarm()].
#' @family Bayesian VAR
#' @export
prior_minnesota <- function(
  lambda = c(0.2, 0.4, 1e-4, 5),
  alpha = c(2, 0.25, 1, 3),
  soc = c(1, 1, 1e-4, 50),
  sur = c(1, 1, 1e-4, 50),
  psi = NULL,
  delta = 1,
  hyper = c("lambda", "soc", "sur"),
  var = 1e7
) {
  fix <- function(x, nm) {
    if (is.null(x)) {
      return(NULL)
    }
    if (length(x) == 1L) {
      x <- c(x, 1, 1e-4, max(50, 100 * x))
    }
    if (length(x) != 4L) {
      rlang::abort(paste0("`", nm, "` must be c(mode, sd, min, max)"))
    }
    stats::setNames(x, c("mode", "sd", "min", "max"))
  }
  if (identical(hyper, "full")) {
    hyper <- c(
      "lambda",
      "alpha",
      if (!is.null(soc)) "soc",
      if (!is.null(sur)) "sur",
      "psi"
    )
  }
  structure(
    list(
      lambda = fix(lambda, "lambda"),
      alpha = fix(alpha, "alpha"),
      soc = fix(soc, "soc"),
      sur = fix(sur, "sur"),
      psi = psi,
      delta = delta,
      hyper = hyper,
      var = var
    ),
    class = "prior_minnesota"
  )
}

gamma_coef <- function(mode, sd) {
  r <- mode^2 / sd^2
  k <- (2 + r + sqrt((4 + r) * r)) / 2
  c(k = k, theta = sqrt(sd^2 / k))
}

#' Bayesian VAR
#'
#' Estimates a Bayesian VAR. The default conjugate Minnesota prior is sampled *exactly*
#' for fixed hyperparameters and with a compiled Metropolis-within-direct-draws sampler for
#' hierarchical hyperparameters (Giannone-Lenza-Primiceri 2015). Everything inside the
#' MCMC loop is C++.
#'
#' @inheritParams varm
#' @param prior [prior_minnesota()] (conjugate, exact/hierarchical) or [prior_triangular()]
#'   (equation-by-equation sampler: Minnesota, SSVS or horseshoe, optional stochastic volatility).
#' @param ndraw number of retained posterior draws.
#' @param nburn burn-in iterations (hierarchical sampler only).
#' @param scale initial scale of the random-walk proposal.
#' @param adapt adapt the proposal scale during burn-in towards 25-45% acceptance.
#' @param sv stochastic volatility (random-walk log-volatility, Kim-Shephard-Chib mixture) --
#'   only with [prior_triangular()].
#' @param thin keep every `thin`-th draw (triangular sampler).
#' @param sv_prior `c(a, b, V0)`: inverse-gamma(`a`, `b`) prior on the log-volatility innovation
#'   variance (default mean 0.01) and prior variance `V0` of the initial log-volatility.
#' @param store_h store the full volatility draws (`ndraw x T x K`); the posterior mean path is always kept.
#' @param vars <[`tidy-select`][tidyselect::language]> endogenous variables
#'   (default: all numeric columns of a data frame).
#' @return a `bvarm` object. Use [irf()], [fevd()], [hd()], [forecast()], [stability()];
#'   `$draws` is a [draws] object, `$hyper` the retained hyperparameter draws, `$logml` the
#'   log marginal likelihood at the mode (or fixed values).
#' @family Bayesian VAR
#' @examples
#' \donttest{
#' fit <- bvarm(yraw, p = 4, ndraw = 2000, nburn = 1000)
#' fit
#' tidy(fit)
#' glance(fit)
#' autoplot(fit, "hyper")
#' plot(irf(fit, horizon = 16))
#' }
#' @export
bvarm <- function(
  Y,
  p = 4L,
  prior = prior_minnesota(),
  ndraw = 5000L,
  nburn = 2000L,
  const = TRUE,
  trend = FALSE,
  exog = NULL,
  scale = 0.5,
  adapt = TRUE,
  sv = FALSE,
  thin = 1L,
  sv_prior = c(a = 5, b = 0.04, V0 = 10),
  store_h = FALSE,
  vars = NULL
) {
  cl <- match.call()
  Y <- prep_data(Y, rlang::enquo(vars))$Y
  if (inherits(prior, "prior_triangular")) {
    return(bvarm_tri(
      Y,
      p,
      prior,
      ndraw,
      nburn,
      thin,
      const,
      trend,
      exog,
      sv,
      sv_prior,
      store_h,
      cl
    ))
  }
  if (sv) {
    rlang::abort("`sv = TRUE` requires `prior = prior_triangular()`")
  }
  if (!inherits(prior, "prior_minnesota")) {
    rlang::abort("`prior` must be prior_minnesota() or prior_triangular()")
  }
  p <- as.integer(p)
  K <- ncol(Y)
  vn <- colnames(Y)
  det <- make_det(nrow(Y), const, trend, exog)
  ols <- cpp_var_ols(Y, p, det)
  m <- ncol(ols$X)
  b <- matrix(0, m, K)
  b[cbind(seq_len(K), seq_len(K))] <- rep_len(prior$delta, K)
  spec <- niw_spec(prior, Y, p, vn)
  args <- list(
    Y = ols$Y,
    X = ols$X,
    p = p,
    d = ncol(det),
    has_const = const,
    b = b,
    vard = prior$var
  )
  fit <- if (any(spec$isfree)) {
    niw_hierarchical(spec, args, ndraw, nburn, scale, adapt)
  } else {
    niw_fixed_draws(spec, args, ndraw)
  }
  dimnames(fit$out$B) <- list(coef_names(vn, p, colnames(det)), vn, NULL)
  draws <- structure(
    list(
      B = fit$out$B,
      Sigma = fit$out$Sigma,
      p = p,
      K = K,
      varnames = vn,
      method = "niw-minnesota"
    ),
    class = "varm_draws"
  )
  new_bvarm(
    draws,
    Y,
    det,
    p,
    K,
    vn,
    list(
      prior = prior,
      hyper = fit$hyper,
      hyper_mode = fit$hmode,
      accept = fit$accept,
      logml = fit$logml,
      call = cl
    )
  )
}

# Hyperparameter vector (start values, free flags, bounds, gamma hyperpriors) of
# the conjugate Minnesota prior; layout lambda, alpha, soc, sur, psi_1..psi_K
niw_spec <- function(prior, Y, p, vn) {
  K <- ncol(Y)
  psi <- prior$psi
  if (is.null(psi)) {
    psi <- ar_resid_var(Y, p)
  }
  hp <- prior[c("lambda", "alpha", "soc", "sur")]
  hyper <- unique(prior$hyper)
  unknown <- setdiff(hyper, c("lambda", "alpha", "soc", "sur", "psi"))
  if (length(unknown)) {
    rlang::abort(paste0(
      "unknown hyperparameter(s): ",
      paste(unknown, collapse = ", ")
    ))
  }
  # a dropped soc/sur dummy is encoded by `off` values in its slot
  pick_or <- function(x, field, off) if (is.null(x)) off else x[field]
  psi_free <- "psi" %in% hyper
  h0 <- c(
    hp$lambda["mode"],
    hp$alpha["mode"],
    pick_or(hp$soc, "mode", -1),
    pick_or(hp$sur, "mode", -1),
    psi
  )
  names(h0) <- c("lambda", "alpha", "soc", "sur", paste0("psi.", vn))
  gamma_prior <- function(x) {
    if (is.null(x)) c(1, 1) else gamma_coef(x["mode"], x["sd"])
  }
  list(
    h0 = h0,
    isfree = c(
      lambda = "lambda" %in% hyper,
      alpha = "alpha" %in% hyper,
      soc = "soc" %in% hyper && !is.null(hp$soc),
      sur = "sur" %in% hyper && !is.null(hp$sur),
      rep(psi_free, K)
    ),
    lo = c(
      hp$lambda["min"],
      hp$alpha["min"],
      pick_or(hp$soc, "min", -Inf),
      pick_or(hp$sur, "min", -Inf),
      if (psi_free) psi / 100 else rep(0, K)
    ),
    hi = c(
      hp$lambda["max"],
      hp$alpha["max"],
      pick_or(hp$soc, "max", Inf),
      pick_or(hp$sur, "max", Inf),
      if (psi_free) psi * 100 else rep(Inf, K)
    ),
    gp = rbind(
      gamma_coef(hp$lambda["mode"], hp$lambda["sd"]),
      gamma_coef(hp$alpha["mode"], hp$alpha["sd"]),
      gamma_prior(hp$soc),
      gamma_prior(hp$sur)
    ),
    ig = c(0.004, 0.004)
  )
}

# Log posterior (or only the log marginal likelihood) of the hyperparameters
niw_logpost <- function(h, spec, args, ml_only = FALSE) {
  do.call(
    cpp_niw_logpost,
    c(
      list(h = h),
      args,
      list(
        isfree = as.integer(spec$isfree),
        gp = spec$gp,
        lo = unname(spec$lo),
        hi = unname(spec$hi),
        ig = spec$ig,
        eval_only_ml = ml_only
      )
    )
  )
}

niw_fixed_draws <- function(spec, args, ndraw) {
  out <- do.call(
    cpp_niw_fixed,
    c(list(h = unname(spec$h0)), args, list(ndraw = as.integer(ndraw)))
  )
  list(out = out, hyper = NULL, accept = NA, logml = out$logml, hmode = spec$h0)
}

# Posterior mode by L-BFGS-B, then the Metropolis-within-direct-draws sampler
niw_hierarchical <- function(spec, args, ndraw, nburn, scale, adapt) {
  h0 <- spec$h0
  fi <- which(spec$isfree)
  objective <- function(x) {
    h <- h0
    h[fi] <- x
    niw_logpost(unname(h), spec, args)
  }
  opt <- stats::optim(
    h0[fi],
    objective,
    method = "L-BFGS-B",
    lower = spec$lo[fi],
    upper = spec$hi[fi],
    hessian = TRUE,
    control = list(fnscale = -1, maxit = 500)
  )
  hmode <- h0
  hmode[fi] <- opt$par
  out <- do.call(
    cpp_niw_hier,
    c(
      list(h0 = unname(hmode)),
      args,
      list(
        isfree = as.integer(spec$isfree),
        gp = spec$gp,
        lo = unname(spec$lo),
        hi = unname(spec$hi),
        ig = spec$ig,
        prop = proposal_cov(opt),
        niter = as.integer(ndraw),
        nburn = as.integer(nburn),
        scale = scale,
        adapt = adapt,
        acc_lo = 0.25,
        acc_hi = 0.45
      )
    )
  )
  hyper <- out$hyper[, fi, drop = FALSE]
  colnames(hyper) <- names(h0)[fi]
  list(
    out = out,
    hyper = hyper,
    accept = out$accept,
    logml = niw_logpost(unname(hmode), spec, args, TRUE),
    hmode = hmode
  )
}

# Random-walk proposal covariance: inverse Hessian at the mode, or a diagonal
# fallback when it is not positive definite
proposal_cov <- function(opt) {
  prop <- tryCatch(solve(-opt$hessian), error = function(e) NULL)
  bad <- is.null(prop) ||
    !all(is.finite(prop)) ||
    any(eigen(prop, symmetric = TRUE)$values <= 0)
  if (bad) {
    prop <- diag(pmax(abs(opt$par), 1e-3)^2 * 0.01, length(opt$par))
  }
  (prop + t(prop)) / 2
}

ar_resid_var <- function(Y, p) {
  vapply(
    seq_len(ncol(Y)),
    function(i) {
      y <- Y[, i]
      n <- length(y)
      X <- cbind(1, stats::embed(y, p + 1)[, -1, drop = FALSE])
      yy <- y[-seq_len(p)]
      r <- yy - X %*% solve(crossprod(X), crossprod(X, yy))
      sum(r^2) / (n - p)
    },
    numeric(1)
  )
}

new_bvarm <- function(draws, Y, det, p, K, vn, extra = list()) {
  structure(
    c(
      list(
        draws = draws,
        Y = Y,
        det = det,
        p = p,
        K = K,
        varnames = vn,
        detnames = colnames(det),
        nobs = nrow(Y) - p
      ),
      extra
    ),
    class = "bvarm"
  )
}

#' @export
print.bvarm <- function(x, digits = 4, ...) {
  N <- dim(x$draws$B)[3]
  cat(
    "<bvarm> ",
    x$draws$method,
    " BVAR(",
    x$p,
    "), K = ",
    x$K,
    ", ",
    N,
    " posterior draws\n",
    sep = ""
  )
  if (!is.null(x$logml)) {
    cat(
      "log marginal likelihood: ",
      format(x$logml, digits = 8),
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$hyper)) {
    cat(
      "hyperparameters (posterior mean [5%, 95%]); MH acceptance ",
      format(x$accept, digits = 3),
      "\n",
      sep = ""
    )
    q <- t(apply(x$hyper, 2, function(z) {
      c(mean(z), stats::quantile(z, c(0.05, 0.95)))
    }))
    colnames(q) <- c("mean", "5%", "95%")
    print(round(q, digits))
  }
  cat("\nposterior mean of B:\n")
  print(round(coef(x), digits))
  invisible(x)
}

#' @export
coef.bvarm <- function(object, ...) apply(object$draws$B, c(1, 2), mean)

#' @export
summary.bvarm <- function(object, ...) {
  B <- object$draws$B
  list(
    mean = apply(B, c(1, 2), mean),
    sd = apply(B, c(1, 2), stats::sd),
    q05 = apply(B, c(1, 2), stats::quantile, 0.05),
    q95 = apply(B, c(1, 2), stats::quantile, 0.95),
    Sigma = apply(object$draws$Sigma, c(1, 2), mean)
  )
}

#' @export
nobs.bvarm <- function(object, ...) object$nobs

#' @export
draws.bvarm <- function(object, ...) object$draws

#' @rdname dynamics
#' @export
irf.bvarm <- function(object, ..., vol = c("last", "mean")) {
  vol <- match.arg(vol)
  d <- object$draws
  if (vol == "mean" && !is.null(object$Sigma_mean)) {
    d$Sigma <- object$Sigma_mean
  }
  irf.varm_draws(d, ...)
}
#' @rdname dynamics
#' @export
fevd.bvarm <- function(object, ..., vol = c("last", "mean")) {
  irf.bvarm(object, ..., vol = vol, .fevd = TRUE)
}
#' @export
stability.bvarm <- function(object) stability.varm_draws(object$draws)

#' @rdname dynamics
#' @export
forecast.bvarm <- function(
  object,
  h = 8L,
  exog_future = NULL,
  level = c(0.68, 0.90),
  ...
) {
  p <- object$p
  Tn <- nrow(object$Y)
  det <- future_det(object, h, exog_future)
  Ylast <- object$Y[(Tn - p + 1):Tn, , drop = FALSE]
  d <- object$draws
  f <- if (isTRUE(object$sv)) {
    cpp_forecast_sv(
      d$B,
      object$A,
      object$h_last,
      object$sigma_h,
      p,
      Ylast,
      det,
      h
    )
  } else {
    cpp_forecast_draws(d$B, d$Sigma, p, Ylast, det, h)
  }
  probs <- level_probs(level)
  bands <- aperm(apply(f, c(2, 3), stats::quantile, probs = probs), c(2, 3, 1))
  dimnames(bands) <- list(NULL, object$varnames, format(probs))
  mean <- apply(f, c(2, 3), mean)
  colnames(mean) <- object$varnames
  structure(
    list(
      mean = mean,
      bands = bands,
      draws = f,
      h = h,
      varnames = object$varnames,
      level = level
    ),
    class = "varm_forecast"
  )
}

#' @rdname dynamics
#' @param ndraw number of posterior draws used for the historical decomposition (median reported).
#' @export
hd.bvarm <- function(object, id = id_chol(), ndraw = 100L, ...) {
  K <- object$K
  p <- object$p
  d <- object$draws
  N <- dim(d$B)[3]
  idx <- unique(round(seq(1, N, length.out = min(ndraw, N))))
  Ye <- object$Y[-seq_len(p), , drop = FALSE]
  X <- cpp_var_ols(object$Y, p, object$det)$X
  acc <- array(0, c(nrow(Ye), K, K, length(idx)))
  for (k in seq_along(idx)) {
    B <- matrix(d$B[,, idx[k]], nrow(d$B), K)
    S <- matrix(d$Sigma[,, idx[k]], K, K)
    P <- structural_irf(
      array(B, c(dim(B), 1)),
      array(S, c(dim(S), 1)),
      K,
      p,
      0L,
      id,
      object$varnames
    )
    P <- matrix(if (id$type == "sign") P$P[,, 1] else P$irf[,, 1, 1], K, K)
    acc[,,, k] <- cpp_hd(B, Ye - X %*% B, P, K, p)
  }
  contrib <- apply(acc, 1:3, stats::median)
  dimnames(contrib) <- list(NULL, object$varnames, object$varnames)
  structure(
    list(contributions = contrib, actual = Ye, varnames = object$varnames),
    class = "varm_hd"
  )
}

#' Triangular-sampler priors (large BVAR, shrinkage, stochastic volatility)
#'
#' Prior for the equation-by-equation sampler of [bvarm()] (Carriero, Clark & Marcellino 2019
#' style): structural equation \eqn{i} regresses \eqn{y_i} on lags, deterministics and the
#' contemporaneous \eqn{y_1,\dots,y_{i-1}}. Variables are therefore *ordered*: the reduced-form
#' fit is order-dependent in the prior (not in the likelihood).
#'
#' The base prior variance of the coefficient on lag \eqn{l} of variable \eqn{j} in equation \eqn{i}
#' is \eqn{\frac{s_i^2}{s_j^2 l^{\alpha}}} (\eqn{s} = AR residual sd), multiplied by
#' * `minnesota`: \eqn{\lambda_1^2} (own lags) or \eqn{\lambda_1^2\lambda_2^2} (other lags);
#' * `ssvs`: slab \eqn{\lambda_1^2} or spike \eqn{\rho^2\lambda_1^2} by a Bernoulli inclusion indicator;
#' * `horseshoe`: \eqn{\tau_i^2\lambda_{ij}^2} with half-Cauchy global/local scales (data-driven shrinkage).
#'
#' @param type `"minnesota"`, `"ssvs"` or `"horseshoe"`.
#' @param lambda1 overall tightness (also the SSVS slab scale).
#' @param lambda2 cross-variable shrinkage factor (Minnesota).
#' @param alpha lag-decay exponent.
#' @param lambda3 prior sd multiplier (times `s_i`) of deterministic coefficients.
#' @param lambda5 prior sd multiplier of the contemporaneous coefficients \eqn{a_{ij}}.
#' @param delta prior mean of the own first lag (1 random walk, 0 stationary).
#' @param rho SSVS spike-to-slab sd ratio.
#' @param inclusion SSVS prior inclusion probability.
#' @param hom prior `c(a0, b0)` for constant variances when `sv = FALSE`: IG(`a0`, `b0 s_i^2`).
#' @return a `prior_triangular` object.
#' @family Bayesian VAR
#' @export
prior_triangular <- function(
  type = c("minnesota", "ssvs", "horseshoe"),
  lambda1 = 0.2,
  lambda2 = 0.5,
  alpha = 2,
  lambda3 = 100,
  lambda5 = 10,
  delta = 1,
  rho = 0.1,
  inclusion = 0.5,
  hom = c(3, 2)
) {
  type <- match.arg(type)
  structure(
    list(
      type = type,
      lambda1 = lambda1,
      lambda2 = lambda2,
      alpha = alpha,
      lambda3 = lambda3,
      lambda5 = lambda5,
      delta = delta,
      rho = rho,
      inclusion = inclusion,
      hom = hom
    ),
    class = "prior_triangular"
  )
}

bvarm_tri <- function(
  Y,
  p,
  prior,
  ndraw,
  nburn,
  thin,
  const,
  trend,
  exog,
  sv,
  sv_prior,
  store_h,
  cl
) {
  Y <- as_data_matrix(Y)
  p <- as.integer(p)
  K <- ncol(Y)
  vn <- colnames(Y)
  det <- make_det(nrow(Y), const, trend, exog)
  ols <- cpp_var_ols(Y, p, det)
  sc <- sqrt(ar_resid_var(Y, p))
  pt <- c(minnesota = 0L, ssvs = 1L, horseshoe = 2L)[[prior$type]]
  pr <- c(
    prior$lambda1,
    prior$lambda2,
    prior$alpha,
    prior$lambda3,
    prior$lambda5,
    prior$rho,
    prior$inclusion
  )
  out <- cpp_bvar_tri(
    ols$Y,
    ols$X,
    p,
    ncol(det),
    pt,
    sv,
    sc,
    rep_len(prior$delta, K),
    pr,
    unname(as.numeric(sv_prior[c(1, 2, 3)])),
    prior$hom,
    as.integer(ndraw),
    as.integer(nburn),
    as.integer(thin),
    store_h
  )
  rn <- coef_names(vn, p, colnames(det))
  dimnames(out$B) <- list(rn, vn, NULL)
  dimnames(out$h_mean) <- list(NULL, vn)
  d <- structure(
    list(
      B = out$B,
      Sigma = out$Sigma,
      p = p,
      K = K,
      varnames = vn,
      method = paste0("triangular-", prior$type, if (sv) "-sv")
    ),
    class = "varm_draws"
  )
  new_bvarm(
    d,
    Y,
    det,
    p,
    K,
    vn,
    list(
      prior = prior,
      sv = sv,
      A = out$A,
      Sigma_mean = out$Sigma_mean,
      h_last = out$h_last,
      sigma_h = out$sigma_h,
      h_mean = out$h_mean,
      h_all = out$h_all,
      inclusion = if (prior$type == "ssvs") out$inclusion,
      tau = if (prior$type == "horseshoe") out$tau,
      call = cl
    )
  )
}

#' Posterior volatility paths
#'
#' Posterior mean of the stochastic volatilities \eqn{\exp(h_{it}/2)} (standard deviation of the
#' structural shock of equation i) from a `bvarm` fitted with `sv = TRUE`.
#' @param object a `bvarm` with stochastic volatility.
#' @return `T_eff x K` matrix.
#' @family Bayesian VAR
#' @export
volatility <- function(object) {
  if (is.null(object$h_mean)) {
    rlang::abort("object has no stochastic volatility")
  }
  exp(object$h_mean / 2)
}
