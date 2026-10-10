#' Local projections
#'
#' Jorda (2005) local projections: for each horizon \eqn{h=0,\dots,H} and each
#' response \eqn{y_k},
#' \deqn{y_{k,t+h} = \alpha_h + \beta_h s_t + \gamma_h' w_t + e_{t+h},}
#' where \eqn{s_t} is the shock, \eqn{w_t} the controls (lags of `Y`, optional
#' contemporaneous variables, extra controls). \eqn{\beta_h} is the impulse
#' response at horizon \eqn{h}. The regressions, 2SLS and HAC standard errors
#' run in C++.
#'
#' @param Y data frame / tibble / matrix (T x K) of responses (and the source of
#'   lag controls). Non-numeric columns (e.g. a date) are dropped.
#' @param shock the shock \eqn{s_t}: a variable of `Y` (tidyselect; observed
#'   shock, unit-impact normalisation) or a numeric vector of length T
#'   (external shock series, e.g. narrative).
#' @param horizon maximum horizon H.
#' @param lags number of lags of `Y` in the controls (and of an external shock).
#' @param contemp <[`tidy-select`][tidyselect::language]> variables of `Y`
#'   entered contemporaneously as controls. For a recursive (Cholesky)
#'   identification with `shock` a variable of `Y`, list the variables ordered
#'   *before* the shock.
#' @param controls optional extra T x q matrix of controls.
#' @param iv optional T x q matrix of instruments for `shock` (LP-IV; Stock &
#'   Watson 2018).
#' @param responses <[`tidy-select`][tidyselect::language]> responses to use
#'   (default all variables of `Y`).
#' @param cumulative <[`tidy-select`][tidyselect::language]> variables whose
#'   response is cumulated (\eqn{\sum_{j=0}^h y_{t+j}}), for variables measured
#'   in differences.
#' @param se `"nw"` Newey-West/Bartlett or `"hc"` Eicker-Huber-White (HC1).
#' @param nw_lag Newey-West truncation lag; default `h + 1`.
#' @param lag_augment Montiel Olea & Plagborg-Moller (2021) lag augmentation:
#'   add one extra lag that is not interpreted and use `"hc"` standard errors,
#'   valid without stationarity assumptions and robust at long horizons.
#' @param level confidence levels for the pointwise bands.
#' @param dof_adjust apply an `n/(n-m)` small-sample factor to the Newey-West
#'   variance.
#' @return an `lpm` object: `$estimate` (`(H+1) x K`), `$se`, `$bands`,
#'   `$nobs`. Use [tidy()][tidy.lpm()] for a tibble and [autoplot()] for a
#'   `ggplot`.
#' @family local projections
#' @examples
#' fit <- lp(yraw, shock = rate, horizon = 12, lags = 2, contemp = c(infl, unemp))
#' fit
#' tidy(fit)
#' @export
lp <- function(
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
  lp_fit(spec)
}

# Shared by lp(), lp_joint(): resolve every argument once -----------------------

lp_spec <- function(
  Y,
  shock,
  contemp,
  responses,
  cumulative,
  horizon,
  lags,
  controls,
  iv,
  se,
  nw_lag,
  lag_augment,
  level,
  dof_adjust,
  call
) {
  Y <- prep_data(Y)$Y
  vn <- colnames(Y)
  if (lag_augment) {
    lags <- lags + 1L
    se <- "hc"
  }
  s <- lp_series(shock, Y, "shock")
  rs <- select_or(responses, vn, "response", seq_along(vn))
  cum <- as.integer(rs %in% select_or(cumulative, vn, "cumulative", integer()))
  d <- lp_design(
    Y,
    s,
    lags,
    select_or(contemp, vn, "contemp", NULL),
    controls
  )
  Z <- if (is.null(iv)) matrix(0, nrow(Y), 0L) else as_data_matrix(iv, "iv")
  H <- as.integer(horizon)
  nwf <- if (is.null(nw_lag)) -1L else as.integer(nw_lag)
  list(
    design = d,
    rs = rs,
    H = H,
    lags = as.integer(lags),
    se = se,
    level = level,
    iv = !is.null(iv),
    call = call,
    engine = list(
      d$Y[, rs, drop = FALSE],
      matrix(d$s, ncol = 1L),
      d$W,
      Z,
      H,
      as.integer(lags),
      cum,
      if (se == "nw") 0L else 1L,
      nwf,
      dof_adjust
    )
  )
}

lp_fit <- function(spec) {
  r <- do.call(cpp_lp, spec$engine)
  H <- spec$H
  vn <- spec$design$vn[spec$rs]
  dn <- list(0:H, vn)
  est <- matrix(r$beta[,, 1], H + 1L, length(vn), dimnames = dn)
  sev <- matrix(r$se[,, 1], H + 1L, length(vn), dimnames = dn)
  probs <- level_probs(spec$level)
  structure(
    list(
      estimate = est,
      se = sev,
      bands = band_array(est, sev, probs),
      probs = probs,
      level = spec$level,
      nobs = r$nobs,
      H = H,
      shock = spec$design$sname,
      responses = vn,
      lags = spec$lags,
      se_type = spec$se,
      iv = spec$iv,
      call = spec$call
    ),
    class = "lpm"
  )
}

# Controls shared by lp(), lp_joint(), lp_state() and lp_smooth(): constant, lags
# of Y (and of an external shock), contemporaneous variables, extra controls.
lp_design <- function(Y, s, lags, contemp = NULL, controls = NULL) {
  if (lags < 1L) {
    rlang::abort("`lags` must be >= 1")
  }
  vn <- colnames(Y)
  W <- cbind(const = 1, lag_matrix(Y, lags, vn))
  if (s$external) {
    W <- cbind(W, lag_matrix(matrix(s$x, ncol = 1L), lags, "shock"))
  }
  if (!is.null(contemp)) {
    W <- cbind(W, Y[, contemp, drop = FALSE])
  }
  if (!is.null(controls)) {
    W <- cbind(W, as_data_matrix(controls, "controls"))
  }
  list(Y = Y, s = s$x, W = W, sname = s$name, vn = vn, n = nrow(Y))
}

# n x (p * ncol) matrix of lags 1..p (NA-padded), columns named `<name>.l<lag>`
lag_matrix <- function(X, p, names = colnames(X)) {
  n <- nrow(X)
  out <- lapply(seq_len(p), function(l) {
    rbind(matrix(NA_real_, l, ncol(X)), X[seq_len(n - l), , drop = FALSE])
  })
  out <- do.call(cbind, out)
  colnames(out) <- as.vector(outer(names, seq_len(p), paste, sep = ".l"))
  out
}

# Data input: drop non-numeric columns, optionally select with tidyselect -------

prep_data <- function(Y, vars = rlang::quo(NULL)) {
  dates <- NULL
  if (is.data.frame(Y)) {
    is_date <- vapply(Y, inherits, NA, what = c("Date", "POSIXt", "yearmon"))
    dates <- if (any(is_date)) Y[[which(is_date)[1]]]
    Y <- Y[vapply(Y, is.numeric, NA)]
  }
  if (!rlang::quo_is_null(vars)) {
    Y <- as.data.frame(Y)[select_vars(vars, names(as.data.frame(Y)), "vars")]
  }
  list(Y = as_data_matrix(Y), dates = dates)
}

# tidyselect with a default for NULL
select_or <- function(quo, names, arg, default) {
  if (rlang::quo_is_null(quo)) {
    return(default)
  }
  select_vars(quo, names, arg)
}

# A series named in Y (tidyselect, one variable) or supplied as a numeric vector
# of length nrow(Y): list(x, name, pos, external)
lp_series <- function(quo, Y, arg) {
  expr <- rlang::quo_get_expr(quo)
  n <- nrow(Y)
  named <- is.character(expr) ||
    (is.symbol(expr) && as.character(expr) %in% colnames(Y))
  val <- if (!named) tryCatch(rlang::eval_tidy(quo), error = function(e) NULL)
  if (is.numeric(val) && length(val) > 1L) {
    if (length(val) != n || anyNA(val)) {
      rlang::abort(
        paste0(
          "external `",
          arg,
          "` must have length nrow(Y), no missing values"
        )
      )
    }
    return(list(x = as.numeric(val), name = arg, pos = NULL, external = TRUE))
  }
  j <- select_vars(quo, colnames(Y), arg)
  if (length(j) != 1L) {
    rlang::abort(paste0("`", arg, "` must select exactly one variable"))
  }
  list(x = Y[, j], name = colnames(Y)[j], pos = j, external = FALSE)
}

# estimate + z * se for every tail probability: array (dim(est), length(probs))
band_array <- function(est, se, probs) {
  z <- stats::qnorm(probs)
  b <- vapply(z, function(zi) est + zi * se, est)
  dimnames(b) <- c(dimnames(est), list(format(probs)))
  b
}

#' @export
print.lpm <- function(x, digits = 4, ...) {
  cat(
    "<lpm> local projection of ",
    paste(x$responses, collapse = ", "),
    " on '",
    x$shock,
    "', horizons 0..",
    x$H,
    " (",
    x$lags,
    " lags, ",
    toupper(x$se_type),
    if (x$iv) ", IV",
    ")\n",
    sep = ""
  )
  rows <- unique(pmin(x$H + 1L, c(1, 2, 5, 9, 13, 21)))
  print(round(x$estimate[rows, , drop = FALSE], digits))
  invisible(x)
}
