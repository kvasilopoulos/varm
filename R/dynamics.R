#' Impulse responses, variance decompositions, historical decompositions, forecasts
#'
#' @param object a `varm`, `bvarm` or `varm_draws` object.
#' @param id identification scheme, see [identification].
#' @param horizon maximum horizon (impact is horizon 0).
#' @param ci for `varm`: `"bootstrap"`, `"wild"`, `"wild_mammen"`, `"mc"` or `"none"`.
#' @param nboot number of bootstrap / Monte Carlo draws.
#' @param level credible/confidence levels (e.g. `c(0.68, 0.90)`); bands are equal-tailed percentiles.
#' @param bias_correct Kilian bias correction (bootstrap only).
#' @param normalise optional `list(shock =, variable =, size =)`: rescale each shock so that the
#'   impact response of `variable` equals `size` (e.g. a 100bp policy-rate shock).
#' @param cumulative character/integer vector of variables whose responses are cumulated over
#'   the horizon (use for variables entered in differences).
#' @param vol for `bvarm` fits with stochastic volatility: use the `"last"` (final-period) or
#'   `"mean"` innovation covariance.
#' @param .fevd internal; compute the variance decomposition instead of the IRF.
#' @param ... passed to methods.
#' @return `varm_irf` / `varm_fevd` objects with `$estimate` (`[K, K, H+1]`, response x shock x horizon;
#'   median across draws when draws exist) and `$bands` (`[K, K, H+1, 2*length(level)]`).
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

level_probs <- function(level) {
  sort(unique(c((1 - level) / 2, 1 - (1 - level) / 2)))
}

cum_h <- function(a, idx) {
  for (i in idx) {
    s <- a[i, , , , drop = FALSE]
    dim(s) <- dim(a)[-1] # shock, h, N
    a[i, , , ] <- aperm(apply(s, c(1, 3), cumsum), c(2, 1, 3))
  }
  a
}

# Fry & Pagan (2011) median-target: single rotation closest to the pointwise median
fry_pagan <- function(a4) {
  d <- dim(a4)
  a <- matrix(a4, ncol = d[4])
  med <- apply(a, 1, stats::median)
  s <- apply(a, 1, stats::sd)
  s[s == 0] <- 1
  which.min(colSums(((a - med) / s)^2))
}

summarise_irf <- function(a4, level, type = "irf") {
  d <- dim(a4)
  K <- d[1]
  H1 <- d[3]
  N <- d[4]
  probs <- level_probs(level)
  q <- cpp_quantiles(array(a4, c(K, K, H1 * N)), H1, c(0.5, probs))
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
  if (!is.null(cumulative)) {
    a4 <- cum_h(a4, resolve_var(cumulative, varnames, "cumulative variable"))
  }
  if (!is.null(normalise)) {
    i <- resolve_var(normalise$variable, varnames)
    j <- resolve_var(normalise$shock, varnames, "shock")
    sz <- if (is.null(normalise$size)) 1 else normalise$size
    a4 <- array(
      cpp_normalise(
        array(a4, c(K, K, dim(a4)[3] * dim(a4)[4])),
        dim(a4)[3],
        i - 1L,
        j - 1L,
        sz
      ),
      dim(a4)
    )
  }
  if (fevd) {
    a4 <- array(
      cpp_fevd_many(array(a4, c(K, K, dim(a4)[3] * dim(a4)[4])), dim(a4)[3]),
      dim(a4)
    )
  }
  s <- if (dim(a4)[4] > 1L) {
    summarise_irf(a4, level)
  } else {
    list(estimate = array(a4[,,, 1], dim(a4)[1:3]), bands = NULL, probs = NULL)
  }
  dn <- list(response = varnames, shock = varnames, h = 0:H)
  dimnames(s$estimate) <- dn
  if (!is.null(s$bands)) {
    dimnames(s$bands) <- c(dn, list(prob = format(s$probs)))
  }
  structure(
    c(
      s,
      list(
        draws = a4,
        K = K,
        H = H,
        id = id,
        varnames = varnames,
        level = level,
        draw = draw
      ),
      extra
    ),
    class = c(if (fevd) "varm_fevd" else "varm_irf", "varm_dyn")
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
  si <- structural_irf(
    object$B,
    object$Sigma,
    object$K,
    object$p,
    horizon,
    id,
    object$varnames
  )
  dyn_object(
    si$irf,
    object$K,
    object$p,
    horizon,
    id,
    object$varnames,
    level,
    normalise,
    cumulative,
    fevd = .fevd,
    draw = si$draw,
    extra = list(accept_rate = si$accept_rate)
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
  d <- if (ci == "none") {
    point_draws(object)
  } else {
    draws(object, ci, nboot, bias_correct)
  }
  res <- irf.varm_draws(
    d,
    id,
    horizon,
    level,
    normalise,
    cumulative,
    .fevd = .fevd
  )
  if (ci != "none") {
    # point estimate from the sample (not the median of draws)
    pt <- point_draws(object)
    if (id$type == "sign") {
      sp <- structural_irf(
        pt$B,
        pt$Sigma,
        object$K,
        object$p,
        horizon,
        id,
        object$varnames
      )
      a <- dyn_object(
        sp$irf,
        object$K,
        object$p,
        horizon,
        id,
        object$varnames,
        level,
        normalise,
        cumulative,
        .fevd
      )
      res$set_estimate <- a$estimate
    } else {
      a <- irf.varm_draws(
        pt,
        id,
        horizon,
        level,
        normalise,
        cumulative,
        .fevd = .fevd
      )
      dn <- dimnames(res$estimate)
      res$estimate <- a$estimate
      dimnames(res$estimate) <- dn
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

#' @rdname dynamics
#' @export
hd.varm <- function(object, id = id_chol(), ...) {
  K <- object$K
  p <- object$p
  pt <- point_draws(object)
  if (id$type == "sign") {
    si <- structural_irf(pt$B, pt$Sigma, K, p, 0L, id, object$varnames)
    P <- matrix(si$P[,, 1], K, K)
  } else {
    P <- matrix(
      structural_irf(pt$B, pt$Sigma, K, p, 0L, id, object$varnames)$irf[,, 1, 1],
      K,
      K
    )
  }
  contrib <- cpp_hd(object$B, object$U, P, K, p)
  dimnames(contrib) <- list(NULL, object$varnames, object$varnames)
  Yd <- object$Y[-seq_len(p), , drop = FALSE]
  structure(
    list(
      contributions = contrib,
      baseline = Yd - apply(contrib, c(1, 2), sum),
      actual = Yd,
      P = P,
      varnames = object$varnames
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
  out <- list(mean = mean, h = h, varnames = object$varnames, level = level)
  if (ci != "none") {
    d <- draws(object, if (ci == "mc") "mc" else ci, nboot)
    f <- cpp_forecast_draws(d$B, d$Sigma, p, Ylast, det, h)
    probs <- level_probs(level)
    out$bands <- aperm(
      apply(f, c(2, 3), stats::quantile, probs = probs),
      c(2, 3, 1)
    )
    dimnames(out$bands) <- list(NULL, object$varnames, format(probs))
    out$draws <- f
  }
  structure(out, class = "varm_forecast")
}

future_det <- function(object, h, exog_future) {
  nd <- ncol(object$det)
  Tn <- nrow(object$Y)
  if (nd == 0L) {
    return(matrix(0, h, 0L))
  }
  cols <- object$detnames
  out <- matrix(0, h, nd, dimnames = list(NULL, cols))
  if ("const" %in% cols) {
    out[, "const"] <- 1
  }
  if ("trend" %in% cols) {
    out[, "trend"] <- Tn + seq_len(h)
  }
  ex <- setdiff(cols, c("const", "trend"))
  if (length(ex)) {
    if (is.null(exog_future)) {
      stop(
        "`exog_future` is required: the model has exogenous regressors",
        call. = FALSE
      )
    }
    exog_future <- as_data_matrix(exog_future, "exog_future")
    if (nrow(exog_future) < h || ncol(exog_future) != length(ex)) {
      stop("`exog_future` must be h x q", call. = FALSE)
    }
    out[, ex] <- exog_future[seq_len(h), ]
  }
  out
}
