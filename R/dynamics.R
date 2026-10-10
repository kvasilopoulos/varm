#' Impulse responses, variance decompositions, historical decompositions, forecasts
#'
#' @param object a `varm`, `bvarm` or `varm_draws` object.
#' @param id identification scheme, see [identification].
#' @param horizon maximum horizon (impact is horizon 0).
#' @param ci for `varm`: `"bootstrap"`, `"wild"`, `"wild_mammen"`, `"mc"` or `"none"`.
#' @param nboot number of bootstrap / Monte Carlo draws.
#' @param level credible/confidence levels (e.g. `c(0.68, 0.90)`); bands are equal-tailed percentiles.
#' @param bias_correct Kilian bias correction (bootstrap only).
#' @param normalise optional [normalise_to()] specification (the older
#'   `list(shock =, variable =, size =)` still works): rescale each shock so that
#'   the impact response of `variable` equals `size` (e.g. a 100bp policy-rate shock).
#' @param cumulative <[`tidy-select`][tidyselect::language]> variables whose
#'   responses are cumulated over the horizon (use for variables entered in
#'   differences), e.g. `cumulative = c(infl, unemp)`.
#' @param vol for `bvarm` fits with stochastic volatility: use the `"last"` (final-period) or
#'   `"mean"` innovation covariance.
#' @param .fevd internal; compute the variance decomposition instead of the IRF.
#' @param ... passed to methods.
#' @return `varm_irf` / `varm_fevd` objects with `$estimate` (`[K, K, H+1]`, response x shock x horizon;
#'   median across draws when draws exist) and `$bands` (`[K, K, H+1, 2*length(level)]`).
#'   Use [tidy()][tidy.varm_irf] for a tibble and [autoplot()][autoplot.varm_irf] for a plot.
#' @family dynamics
#' @examples
#' fit <- varm(yraw, p = 2)
#' r <- irf(fit, horizon = 8, ci = "mc", nboot = 100, cumulative = infl)
#' r
#' tidy(r)
#' ggplot2::autoplot(r)
#'
#' # a 100bp policy-rate shock
#' irf(fit, horizon = 4, ci = "none", normalise = normalise_to(rate, rate, 1))
#'
#' fevd(fit, horizon = 8, ci = "none")
#' forecast(fit, h = 4)
#' @name dynamics
NULL

#' @rdname dynamics
#' @export
irf <- function(object, ...) UseMethod("irf")
#' @rdname dynamics
#' @export
fevd <- function(object, ...) UseMethod("fevd")
#' @rdname dynamics
#' @export
hd <- function(object, ...) UseMethod("hd")
#' @rdname dynamics
#' @export
forecast <- function(object, ...) UseMethod("forecast")

#' Normalise the impact of a shock
#'
#' Rescales each structural shock so that the impact response of `variable` to
#' `shock` equals `size`. Pass the result to the `normalise` argument of [irf()].
#'
#' @param shock,variable <[`tidy-select`][tidyselect::language]> a single
#'   variable each: the shock being rescaled and the variable whose impact
#'   response is fixed.
#' @param size impact response of `variable` (default 1).
#' @return a `varm_normalise` specification.
#' @family dynamics
#' @examples
#' normalise_to(rate, rate, size = 2)
#' @export
normalise_to <- function(shock, variable, size = 1) {
  structure(
    list(
      shock = rlang::enquo(shock),
      variable = rlang::enquo(variable),
      size = size
    ),
    class = "varm_normalise"
  )
}

# helpers on IRF arrays [K, K, H+1, N] ----------------------------------------

# flatten the horizon and draw dimensions for the compiled routines
flat3 <- function(a) array(a, c(dim(a)[1:2], prod(dim(a)[3:4])))

# cumulate over the horizon the responses of the variables `idx`
cum_h <- function(a, idx) {
  for (i in idx) {
    s <- apply(a[i, , , , drop = FALSE], c(2, 4), cumsum)
    a[i, , , ] <- aperm(s, c(2, 1, 3))
  }
  a
}

normalise_irf <- function(a4, spec, varnames) {
  i <- select_one(spec$variable, varnames, "variable")
  j <- select_one(spec$shock, varnames, "shock")
  size <- spec$size %||% 1
  H1 <- dim(a4)[3]
  array(cpp_normalise(flat3(a4), H1, i - 1L, j - 1L, size), dim(a4))
}

# Fry & Pagan (2011) median-target: single rotation closest to the pointwise median
fry_pagan <- function(a4) {
  a <- matrix(a4, ncol = dim(a4)[4])
  med <- apply(a, 1, stats::median)
  s <- apply(a, 1, stats::sd)
  s[s == 0] <- 1
  which.min(colSums(((a - med) / s)^2))
}

# pointwise median and percentile bands across draws
summarise_irf <- function(a4, level) {
  K <- dim(a4)[1]
  H1 <- dim(a4)[3]
  if (dim(a4)[4] == 1L) {
    return(list(estimate = array(a4, c(K, K, H1)), bands = NULL, probs = NULL))
  }
  probs <- level_probs(level)
  q <- cpp_quantiles(flat3(a4), H1, c(0.5, probs))
  q <- array(q, c(K, K, H1, length(probs) + 1L))
  list(
    estimate = array(q[,,, 1], c(K, K, H1)),
    bands = q[,,, -1, drop = FALSE],
    probs = probs
  )
}

dyn_object <- function(
  a4,
  K,
  p,
  H,
  id,
  varnames,
  level,
  normalise,
  cumulative,
  fevd = FALSE,
  draw = NULL,
  extra = list()
) {
  cum <- select_vars(cumulative, varnames, "cumulative variable")
  if (length(cum)) {
    a4 <- cum_h(a4, cum)
  }
  if (!is.null(normalise)) {
    a4 <- normalise_irf(a4, normalise, varnames)
  }
  if (fevd) {
    a4 <- array(cpp_fevd_many(flat3(a4), dim(a4)[3]), dim(a4))
  }
  s <- summarise_irf(a4, level)
  dn <- list(response = varnames, shock = varnames, h = 0:H)
  dimnames(s$estimate) <- dn
  if (!is.null(s$bands)) {
    dimnames(s$bands) <- c(dn, list(prob = format(s$probs)))
  }
  info <- list(
    draws = a4,
    K = K,
    H = H,
    id = id,
    varnames = varnames,
    level = level,
    draw = draw
  )
  structure(
    c(s, info, extra),
    class = c(if (fevd) "varm_fevd" else "varm_irf", "varm_dyn")
  )
}

# structural IRF of a `varm_draws` object (shared by irf() and fevd())
irf_core <- function(
  d,
  id,
  horizon,
  level,
  normalise,
  cumulative,
  fevd = FALSE
) {
  si <- structural_irf(d$B, d$Sigma, d$K, d$p, horizon, id, d$varnames)
  dyn_object(
    si$irf,
    d$K,
    d$p,
    horizon,
    id,
    d$varnames,
    level,
    normalise,
    cumulative,
    fevd = fevd,
    draw = si$draw,
    extra = list(accept_rate = si$accept_rate)
  )
}

#' @rdname dynamics
#' @export
irf.varm_draws <- function(
  object,
  id = id_chol(),
  horizon = 20L,
  level = c(0.68, 0.90),
  normalise = NULL,
  cumulative = NULL,
  ...,
  .fevd = FALSE
) {
  irf_core(
    object,
    id,
    horizon,
    level,
    normalise,
    rlang::enquo(cumulative),
    .fevd
  )
}

#' @rdname dynamics
#' @export
fevd.varm_draws <- function(
  object,
  id = id_chol(),
  horizon = 20L,
  level = c(0.68, 0.90),
  ...
) {
  irf.varm_draws(object, id, horizon, level, ..., .fevd = TRUE)
}

# a single "draw": the point estimate
point_draws <- function(object) {
  structure(
    list(
      B = array(object$B, c(dim(object$B), 1L)),
      Sigma = array(object$Sigma, c(dim(object$Sigma), 1L)),
      p = object$p,
      K = object$K,
      varnames = object$varnames,
      method = "point",
      fit = object
    ),
    class = "varm_draws"
  )
}

#' @rdname dynamics
#' @export
irf.varm <- function(
  object,
  id = id_chol(),
  horizon = 20L,
  ci = c("bootstrap", "wild", "wild_mammen", "mc", "none"),
  nboot = 1000L,
  level = c(0.68, 0.90),
  bias_correct = FALSE,
  normalise = NULL,
  cumulative = NULL,
  ...,
  .fevd = FALSE
) {
  ci <- match.arg(ci)
  cumulative <- rlang::enquo(cumulative)
  d <- if (ci == "none") {
    point_draws(object)
  } else {
    draws(object, ci, nboot, bias_correct)
  }
  res <- irf_core(d, id, horizon, level, normalise, cumulative, .fevd)
  if (ci != "none") {
    # point estimate from the sample (not the median of draws)
    pt <- irf_core(
      point_draws(object),
      id,
      horizon,
      level,
      normalise,
      cumulative,
      .fevd
    )
    if (id$type == "sign") {
      res$set_estimate <- pt$estimate
    } else {
      res$estimate[] <- pt$estimate
    }
  }
  res$method <- ci
  res
}

#' @rdname dynamics
#' @export
fevd.varm <- function(object, id = id_chol(), horizon = 20L, ...) {
  irf.varm(object, id, horizon, ..., .fevd = TRUE)
}

# impact matrix P (response x shock) of the point estimate
impact_matrix <- function(object, id) {
  pt <- point_draws(object)
  K <- object$K
  si <- structural_irf(pt$B, pt$Sigma, K, object$p, 0L, id, object$varnames)
  matrix(if (id$type == "sign") si$P[,, 1] else si$irf[,, 1, 1], K, K)
}

#' @rdname dynamics
#' @export
hd.varm <- function(object, id = id_chol(), ...) {
  K <- object$K
  p <- object$p
  P <- impact_matrix(object, id)
  contrib <- cpp_hd(object$B, object$U, P, K, p)
  dimnames(contrib) <- list(NULL, object$varnames, object$varnames)
  Yd <- object$Y[-seq_len(p), , drop = FALSE]
  structure(
    list(
      contributions = contrib,
      baseline = Yd - apply(contrib, c(1, 2), sum),
      actual = Yd,
      P = P,
      varnames = object$varnames,
      dates = object$dates[-seq_len(p)]
    ),
    class = "varm_hd"
  )
}

#' @rdname dynamics
#' @param h forecast horizon (steps).
#' @param exog_future `h x q` matrix of future exogenous regressors (if `exog` was used).
#' @export
forecast.varm <- function(
  object,
  h = 8L,
  exog_future = NULL,
  ci = c("mc", "bootstrap", "wild", "none"),
  nboot = 1000L,
  level = c(0.68, 0.90),
  ...
) {
  ci <- match.arg(ci)
  K <- object$K
  p <- object$p
  Tn <- nrow(object$Y)
  det <- future_det(object, h, exog_future)
  Ylast <- object$Y[(Tn - p + 1):Tn, , drop = FALSE]
  mean <- cpp_path(object$B, K, p, Ylast, det, matrix(0, h, K))
  colnames(mean) <- object$varnames
  out <- list(
    mean = mean,
    h = h,
    varnames = object$varnames,
    level = level,
    history = object$Y,
    dates = object$dates
  )
  if (ci != "none") {
    d <- draws(object, ci, nboot)
    out$draws <- cpp_forecast_draws(d$B, d$Sigma, p, Ylast, det, h)
    out$bands <- forecast_bands(
      out$draws,
      level_probs(level),
      object$varnames
    )
  }
  structure(out, class = "varm_forecast")
}

future_det <- function(object, h, exog_future) {
  if (ncol(object$det) == 0L) {
    return(matrix(0, h, 0L))
  }
  cols <- object$detnames
  out <- matrix(0, h, length(cols), dimnames = list(NULL, cols))
  if ("const" %in% cols) {
    out[, "const"] <- 1
  }
  if ("trend" %in% cols) {
    out[, "trend"] <- nrow(object$Y) + seq_len(h)
  }
  ex <- setdiff(cols, c("const", "trend"))
  if (length(ex)) {
    out[, ex] <- future_exog(exog_future, h, length(ex))
  }
  out
}

future_exog <- function(exog_future, h, q) {
  if (is.null(exog_future)) {
    rlang::abort(
      "`exog_future` is required: the model has exogenous regressors"
    )
  }
  exog_future <- as_data_matrix(exog_future, "exog_future")
  if (nrow(exog_future) < h || ncol(exog_future) != q) {
    rlang::abort("`exog_future` must be h x q")
  }
  exog_future[seq_len(h), ]
}

# printing ------------------------------------------------------------------

#' @export
print.varm_dyn <- function(x, ...) {
  what <- if (inherits(x, "varm_fevd")) "FEVD" else "IRF"
  bands <- if (!is.null(x$bands)) {
    paste0(", bands: ", paste0(x$level * 100, "%", collapse = "/"))
  }
  cat(
    "<varm ",
    what,
    "> K = ",
    x$K,
    ", horizon 0..",
    x$H,
    ", id = ",
    x$id$type,
    bands,
    "\n",
    sep = ""
  )
  if (!is.na(x$accept_rate %||% NA)) {
    cat(
      "sign-restriction acceptance rate: ",
      format(x$accept_rate, digits = 3),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' @export
print.varm_forecast <- function(x, ...) {
  cat("<varm forecast> h = ", x$h, "\n", sep = "")
  print(round(x$mean, 4))
  invisible(x)
}
