#' Smooth local projections (Barnichon & Brownlees)
#'
#' Imposes that the impulse response is a smooth function of the horizon,
#' \eqn{\beta(h) = \sum_{j=1}^J \theta_j B_j(h)} with a B-spline basis \eqn{B_j} (cubic by
#' default, equally spaced knots over \eqn{0..H}), and estimates all horizons jointly by penalised least squares:
#' \deqn{\hat\theta = \arg\min_\theta \sum_{h=0}^H\sum_t \Big(\tilde y_{t+h}-\sum_j\theta_jB_j(h)\tilde s_{t,h}\Big)^2 +
#'   \lambda\,\theta'D_r'D_r\theta,}
#' where tilde denotes residuals from the horizon-specific regression on the controls (Frisch-Waugh),
#' \eqn{D_r} is the \eqn{r}-th order difference matrix (P-spline penalty, Eilers & Marx 1996) which shrinks
#' \eqn{\beta(h)} towards a polynomial of degree \eqn{r-1} as \eqn{\lambda\to\infty}. \eqn{\lambda}
#' is chosen by K-fold cross-validation over contiguous blocks of time (preserving serial
#' dependence). Standard errors use the HAC sandwich
#' \eqn{(A+\lambda P)^{-1}\hat\Omega(A+\lambda P)^{-1}}, \eqn{\hat\Omega} Newey-West on the scores
#' summed over horizons at each \eqn{t}; they treat the smoothing bias as negligible.
#'
#' **Economic intuition.** Raw LP estimates are noisy and jagged because every horizon is
#' estimated on its own; theory says responses are smooth (hump-shaped, monotone decay). Smoothing
#' borrows strength across neighbouring horizons, trading a small bias (a hump is flattened if
#' \eqn{\lambda} is too large) for a large variance reduction, which is what the LP-vs-VAR
#' bias-variance trade-off is about: smooth LP sits between unrestricted LP and a parametric VAR.
#' *Pitfalls:* sharp features (e.g. an impact jump at \eqn{h=0}, a sudden reversal) are smoothed
#' away; use enough basis functions; the bands are conditional on \eqn{\lambda}, ignoring the CV
#' step; one response at a time.
#'
#' @inheritParams lp
#' @param response <[`tidy-select`][tidyselect::language]> the single response
#'   variable of `Y`.
#' @param cumulative cumulate the response (for variables in differences).
#' @param nbasis number of B-spline basis functions \eqn{J} (default `min(H + 1, 8)`).
#' @param degree spline degree.
#' @param penalty_order order \eqn{r} of the difference penalty.
#' @param lambda a fixed (absolute) penalty; `NULL` = choose by cross-validation over a grid
#'   scaled by \eqn{\mathrm{tr}(A)/J}.
#' @param nfolds number of CV blocks.
#' @param nw_lag Newey-West truncation for the score (default `H + 1`).
#' @return an `lp_smooth` object: `$estimate`/`$se` (`(H+1) x 1`), `$bands`, `$vcov` (compatible
#'   with [lp_supt()]), `$lp`/`$lp_se` (unsmoothed LP), `$lambda`, `$cv` (grid and CV errors), `$edf`.
#' @references Barnichon, R. and Brownlees, C. (2019) Impulse response estimation by smooth local
#'   projections. *Review of Economics and Statistics*.
#' @family local projections
#' @examples
#' fit <- lp_smooth(
#'   yraw,
#'   shock = rate,
#'   response = infl,
#'   horizon = 12,
#'   lags = 2,
#'   contemp = c(infl, unemp)
#' )
#' round(cbind(smooth = fit$estimate[, 1], raw = fit$lp), 3)
#' @export
lp_smooth <- function(
  Y,
  shock,
  response = 1L,
  horizon = 20L,
  lags = 4L,
  contemp = NULL,
  controls = NULL,
  cumulative = FALSE,
  nbasis = NULL,
  degree = 3L,
  penalty_order = 2L,
  lambda = NULL,
  nfolds = 5L,
  nw_lag = NULL,
  level = c(0.68, 0.90)
) {
  H <- as.integer(horizon)
  lags <- as.integer(lags)
  Y <- prep_data(Y)$Y
  vn <- colnames(Y)
  s <- lp_series(rlang::enquo(shock), Y, "shock")
  rj <- select_vars(rlang::enquo(response), vn, "response")
  if (length(rj) != 1L) {
    rlang::abort("`response` must select exactly one variable")
  }
  d <- lp_design(
    Y,
    s,
    lags,
    select_or(rlang::enquo(contemp), vn, "contemp", NULL),
    controls
  )
  J <- if (is.null(nbasis)) min(H + 1L, 8L) else as.integer(nbasis)
  if (J > H + 1L) {
    rlang::abort("`nbasis` cannot exceed horizon + 1")
  }
  degree <- min(as.integer(degree), J - 1L)
  penalty_order <- min(as.integer(penalty_order), J - 1L)
  Bm <- smooth_basis(H, J, degree)
  st <- smooth_stack(d, Y[, rj], H, lags, cumulative)
  X <- Bm[st$h + 1L, , drop = FALSE] * st$s
  P <- crossprod(diff(diag(J), differences = penalty_order))
  cv <- smooth_cv(X, st$y, st$t, d$n - lags, nfolds, P, lambda)
  Mi <- solve(cv$At + cv$lambda * P)
  th <- drop(Mi %*% cv$bt)
  Sc <- rowsum(X * (st$y - drop(X %*% th)), st$t)
  V <- Mi %*%
    nw_lrv(Sc, if (is.null(nw_lag)) H + 1L else as.integer(nw_lag)) %*%
    Mi
  est <- drop(Bm %*% th)
  Vb <- Bm %*% V %*% t(Bm)
  sev <- sqrt(pmax(diag(Vb), 0))
  dn <- list(0:H, vn[rj])
  est <- matrix(est, dimnames = dn)
  sev <- matrix(sev, dimnames = dn)
  probs <- level_probs(level)
  dimnames(Vb) <- list(0:H, 0:H)
  structure(
    list(
      estimate = est,
      se = sev,
      bands = band_array(est, sev, probs),
      probs = probs,
      level = level,
      vcov = stats::setNames(list(Vb), vn[rj]),
      lp = st$lp,
      lp_se = st$lp_se,
      lambda = cv$lambda,
      cv = data.frame(lambda = cv$grid, cv = cv$cv),
      theta = th,
      basis = Bm,
      edf = sum(diag(Mi %*% cv$At)),
      H = H,
      shock = d$sname,
      responses = vn[rj],
      nobs = d$n - lags - 0:H
    ),
    class = c("lp_smooth", "lpm")
  )
}

# B-spline basis (H+1 x J), equally spaced knots over 0..H
smooth_basis <- function(H, J, degree) {
  dx <- H / (J - degree)
  knots <- seq(-degree * dx, H + degree * dx, length.out = J + degree + 1L)
  splines::splineDesign(knots, 0:H, ord = degree + 1L, outer.ok = TRUE)
}

# Frisch-Waugh residuals of the response and the shock at every horizon, stacked,
# with the unsmoothed LP estimate and HAC standard error of each horizon
smooth_stack <- function(d, y, H, lags, cumulative) {
  cs <- cumsum(y)
  one <- function(h) {
    ix <- (lags + 1L):(d$n - h)
    yh <- if (cumulative) cs[ix + h] - cs[ix - 1L] else y[ix + h]
    qw <- qr(d$W[ix, , drop = FALSE])
    yr <- qr.resid(qw, yh)
    sr <- qr.resid(qw, d$s[ix])
    sx <- sum(sr^2)
    b <- sum(sr * yr) / sx
    g <- sr * (yr - b * sr) / sx
    list(
      y = yr,
      s = sr,
      t = ix - lags,
      h = rep(h, length(ix)),
      b = b,
      se = sqrt(nw_lrv(g, h + 1L))
    )
  }
  res <- lapply(0:H, one)
  pull <- function(nm) unlist(lapply(res, `[[`, nm))
  list(
    y = pull("y"),
    s = pull("s"),
    t = pull("t"),
    h = pull("h"),
    lp = pull("b"),
    lp_se = pull("se")
  )
}

# Block K-fold CV over contiguous time blocks; returns the chosen penalty
smooth_cv <- function(X, y, tt, n0, nfolds, P, lambda) {
  fold <- pmin(nfolds, 1L + ((tt - 1L) * nfolds) %/% n0)
  Af <- lapply(seq_len(nfolds), function(f) {
    crossprod(X[fold == f, , drop = FALSE])
  })
  bf <- lapply(seq_len(nfolds), function(f) {
    crossprod(X[fold == f, , drop = FALSE], y[fold == f])
  })
  cf <- vapply(seq_len(nfolds), function(f) sum(y[fold == f]^2), 0)
  At <- Reduce(`+`, Af)
  bt <- Reduce(`+`, bf)
  grid <- if (is.null(lambda)) {
    mean(diag(At)) * 10^seq(-4, 4, length.out = 21)
  } else {
    lambda
  }
  fold_err <- function(f, l) {
    th <- solve(At - Af[[f]] + l * P, bt - bf[[f]])
    cf[f] - 2 * sum(th * bf[[f]]) + drop(crossprod(th, Af[[f]] %*% th))
  }
  cv <- vapply(
    grid,
    function(l) sum(vapply(seq_len(nfolds), fold_err, 0, l = l)),
    0
  )
  best <- if (length(grid) == 1L) 1L else which.min(cv)
  list(lambda = grid[best], grid = grid, cv = cv, At = At, bt = bt)
}
