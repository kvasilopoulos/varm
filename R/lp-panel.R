#' Panel local projections with fixed effects and clustered errors
#'
#' For a panel \eqn{(i,t)}, estimates for each horizon
#' \deqn{y_{i,t+h} = \alpha_{i,h} + \delta_{t,h} + \beta_h s_{it} + \sum_{l=1}^{L}(\phi_{l,h} y_{i,t-l} + \psi_{l,h}s_{i,t-l}) +
#'   \kappa_h' x_{it} + e_{i,t+h}}
#' by within (demeaning) estimation with unit effects \eqn{\alpha_i}, time effects \eqn{\delta_t},
#' both, or neither. With two-way effects on an unbalanced horizon sample the demeaning is iterated
#' (alternating projections), a single pass for balanced panels, so no dummy matrix is ever built.
#' Standard errors cluster on the unit (or time) with the Stata-type small-sample correction
#' \eqn{\frac{G}{G-1}\frac{n-1}{n-k}}, \eqn{k} counting the regressors and the absorbed effects (identical to
#' `lm` with dummies and hand-coded cluster scores). The joint cross-horizon cluster covariance
#' (scores of a cluster summed over time, stacked over horizons) is returned for [lp_supt()].
#'
#' **Economic intuition.** Panels (countries, states, firms) give identification and power: unit
#' effects remove time-invariant heterogeneity and time effects remove common shocks, so the
#' response is identified from cross-sectional differences in exposure to \eqn{s}. Clustering by unit
#' allows arbitrary serial correlation within the unit, which is essential because the LP
#' residual is an MA(\eqn{h}) process. *Pitfalls:* fixed effects with lagged dependent variables
#' induce Nickell bias when \eqn{T} is small (use long panels or drop `y` lags); time effects absorb the
#' average response, so the estimate is a *relative* effect; with few clusters (< 30) the
#' cluster-robust bands are too narrow; heterogeneous responses are averaged with variance weights.
#'
#' @param data a data frame in long format.
#' @param y,shock,unit,time <[`tidy-select`][tidyselect::language]> the columns of
#'   the outcome, the shock, the unit id and the time index (regularly spaced).
#' @param horizon maximum horizon.
#' @param lags number of lags of `y` and `shock` among the controls.
#' @param controls <[`tidy-select`][tidyselect::language]> further
#'   (contemporaneous) control columns.
#' @param fe `"twoway"`, `"unit"`, `"time"` or `"none"` (common intercept).
#' @param cluster `"unit"` or `"time"`.
#' @param cumulative cumulate the outcome over the horizon.
#' @param level confidence levels.
#' @return an `lp_panel` object with `$estimate`, `$se`, `$bands` (`(H+1) x 1`), `$vcov`, `$nobs`, `$nclusters`.
#' @family local projections
#' @examples
#' set.seed(1)
#' d <- expand.grid(id = 1:20, t = 1:40)
#' d$s <- rnorm(nrow(d))
#' d$y <- 0.5 * d$s + rnorm(nrow(d))
#' lp_panel(d, y, s, id, t, horizon = 3, lags = 1)$estimate
#' @export
lp_panel <- function(
  data,
  y,
  shock,
  unit,
  time,
  horizon = 8L,
  lags = 1L,
  controls = NULL,
  fe = c("twoway", "unit", "time", "none"),
  cluster = c("unit", "time"),
  cumulative = FALSE,
  level = c(0.68, 0.90)
) {
  fe <- match.arg(fe)
  cluster <- match.arg(cluster)
  d <- as.data.frame(data)
  nms <- names(d)
  one <- function(quo, arg) {
    j <- select_vars(quo, nms, arg)
    if (length(j) != 1L) {
      rlang::abort(paste0("`", arg, "` must select exactly one column"))
    }
    nms[j]
  }
  y <- one(rlang::enquo(y), "y")
  shock <- one(rlang::enquo(shock), "shock")
  unit <- one(rlang::enquo(unit), "unit")
  time <- one(rlang::enquo(time), "time")
  ctrl <- nms[select_or(rlang::enquo(controls), nms, "controls", integer())]

  H <- as.integer(horizon)
  ix <- panel_index(d[[unit]], d[[time]], H, as.integer(lags))
  yv <- d[[y]]
  X <- panel_regressors(yv, d[[shock]], ix$pos, as.integer(lags), d[ctrl])
  cl <- if (cluster == "unit") ix$u else ix$tm

  est <- se <- nob <- ncl <- numeric(H + 1L)
  Sc <- matrix(0, nrow(d), H + 1L)
  yh <- yv
  for (h in 0:H) {
    if (h > 0) {
      yh <- if (cumulative) yh + yv[ix$pos(h)] else yv[ix$pos(h)]
    }
    f <- panel_horizon(yh, X, fe, ix$u, ix$tm, cl)
    est[h + 1L] <- f$est
    se[h + 1L] <- f$se
    nob[h + 1L] <- f$n
    ncl[h + 1L] <- f$G
    Sc[f$ok, h + 1L] <- f$q
  }
  V <- crossprod(rowsum(Sc, cl))
  dimnames(V) <- list(0:H, 0:H)
  dn <- list(0:H, y)
  est <- matrix(est, dimnames = dn)
  se <- matrix(se, dimnames = dn)
  probs <- level_probs(level)
  structure(
    list(
      estimate = est,
      se = se,
      bands = band_array(est, se, probs),
      probs = probs,
      level = level,
      vcov = stats::setNames(list(V), y),
      nobs = nob,
      nclusters = ncl,
      H = H,
      shock = shock,
      responses = y,
      fe = fe,
      cluster = cluster
    ),
    class = c("lp_panel", "lpm")
  )
}

# Integer unit and period indices and a lead/lag lookup `pos(k)` (NA if absent)
panel_index <- function(unit, time, H, L) {
  u <- as.integer(factor(unit))
  ut <- sort(unique(time))
  step <- if (length(ut) > 1L) min(diff(ut)) else 1
  tm <- as.integer(round((time - min(ut)) / step)) + 1L
  key <- (u - 1L) * (max(tm) + H + L + 5L) + tm
  if (anyDuplicated(key)) {
    rlang::abort("duplicated (unit, time) pairs")
  }
  list(u = u, tm = tm, pos = function(k) match(key + k, key))
}

# Contemporaneous shock, lags of y and shock, extra controls
panel_regressors <- function(y, s, pos, L, extra) {
  X <- cbind(shock = s)
  for (l in seq_len(L)) {
    X <- cbind(X, y[pos(-l)], s[pos(-l)])
  }
  if (L > 0) {
    colnames(X) <- c(
      "shock",
      paste0(c("y.l", "shock.l"), rep(seq_len(L), each = 2L))
    )
  }
  if (ncol(extra) > 0L) {
    X <- cbind(X, as.matrix(extra))
  }
  X
}

# Within estimate of the shock coefficient at one horizon, with cluster scores
panel_horizon <- function(yh, X, fe, u, tm, cl) {
  ok <- stats::complete.cases(yh, X)
  n <- sum(ok)
  g <- switch(
    fe,
    twoway = list(u[ok], tm[ok]),
    unit = list(u[ok]),
    time = list(tm[ok]),
    none = list(rep(1L, n))
  )
  g <- lapply(g, function(v) match(v, unique(v)))
  Mm <- demean_within(cbind(yh[ok], X[ok, , drop = FALSE]), g)
  yd <- Mm[, 1]
  Xd <- Mm[, -1, drop = FALSE]
  Ri <- solve(crossprod(Xd))
  b <- drop(Ri %*% crossprod(Xd, yd))
  e <- yd - drop(Xd %*% b)
  n_u <- length(unique(u[ok]))
  n_t <- length(unique(tm[ok]))
  kfe <- switch(fe, twoway = n_u + n_t - 1L, unit = n_u, time = n_t, none = 1L)
  G <- length(unique(cl[ok]))
  cadj <- G / (G - 1) * (n - 1) / (n - ncol(X) - kfe)
  q <- drop(Xd %*% Ri[, 1]) * e * sqrt(cadj)
  list(
    est = b[1],
    se = sqrt(sum(rowsum(q, cl[ok])^2)),
    n = n,
    G = G,
    ok = ok,
    q = q
  )
}

# iterated within transformation (alternating projections); one pass if balanced
demean_within <- function(M, groups, tol = 1e-10, maxit = 200L) {
  M <- as.matrix(M)
  gs <- lapply(groups, function(g) list(g = g, cnt = tabulate(g)))
  for (it in seq_len(maxit)) {
    mx <- 0
    for (f in gs) {
      m <- rowsum(M, f$g, reorder = TRUE) / f$cnt
      M <- M - m[f$g, , drop = FALSE]
      mx <- max(mx, abs(m))
    }
    if (length(gs) == 1L || mx < tol) {
      break
    }
  }
  M
}
