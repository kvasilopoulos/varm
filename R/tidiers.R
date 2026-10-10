# Tidiers (broom conventions): every method returns a tibble ------------------

#' Tidy a VAR fit
#'
#' Coefficient table of a [varm()] fit in the broom convention, a one-row model
#' summary, and the data augmented with fitted values and residuals.
#'
#' @param x a `varm` fit.
#' @param conf.int include `conf.low` / `conf.high` (t-based).
#' @param conf.level confidence level for `conf.int`.
#' @param data data to augment; defaults to the data used for the fit.
#' @param ... unused.
#' @return `tidy()`: a tibble with `response`, `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value` (and `conf.low`, `conf.high`). `glance()`: one row
#'   with `K`, `p`, `nobs`, `logLik`, `AIC`, `BIC`, `HQ`, `max_root`.
#'   `augment()`: the data plus `.fitted_<var>` and `.resid_<var>` columns (`NA`
#'   for the first `p` rows).
#' @family tidiers
#' @examples
#' fit <- varm(yraw, p = 2)
#' tidy(fit)
#' tidy(fit, conf.int = TRUE)
#' glance(fit)
#' augment(fit)
#' @name tidy.varm
NULL

#' @rdname tidy.varm
#' @export
tidy.varm <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  s <- summary.varm(x)
  out <- tibble::tibble(
    response = rep(x$varnames, each = nrow(x$B)),
    term = rep(rownames(x$B), times = x$K),
    estimate = as.vector(x$B),
    std.error = as.vector(s$se),
    statistic = as.vector(s$tvalue),
    p.value = as.vector(s$pvalue)
  )
  if (conf.int) {
    crit <- stats::qt(1 - (1 - conf.level) / 2, x$nobs - x$npar)
    out$conf.low <- out$estimate - crit * out$std.error
    out$conf.high <- out$estimate + crit * out$std.error
  }
  out
}

#' @rdname tidy.varm
#' @export
glance.varm <- function(x, ...) {
  ll <- logLik(x)
  df <- attr(ll, "df")
  tibble::tibble(
    K = x$K,
    p = x$p,
    nobs = x$nobs,
    logLik = as.numeric(ll),
    AIC = stats::AIC(ll),
    BIC = stats::BIC(ll),
    HQ = -2 * as.numeric(ll) + 2 * df * log(log(x$nobs)),
    max_root = max(stability(x))
  )
}

#' @rdname tidy.varm
#' @export
augment.varm <- function(x, data = NULL, ...) {
  if (is.null(data)) {
    data <- tibble::as_tibble(x$Y)
    if (!is.null(x$dates)) {
      data <- tibble::add_column(data, date = x$dates, .before = 1L)
    }
  }
  if (nrow(data) != nrow(x$Y)) {
    rlang::abort("`data` must have the same number of rows as the fitted data")
  }
  pad <- matrix(NA_real_, x$p, x$K)
  fit <- rbind(pad, fitted(x))
  res <- rbind(pad, resid(x))
  out <- tibble::as_tibble(data)
  for (k in seq_len(x$K)) {
    out[[paste0(".fitted_", x$varnames[k])]] <- fit[, k]
    out[[paste0(".resid_", x$varnames[k])]] <- res[, k]
  }
  out
}

# impulse responses and variance decompositions -------------------------------

#' Tidy impulse responses and variance decompositions
#'
#' @param x a `varm_irf` or `varm_fevd` object.
#' @param bands include the percentile bands (one row per cell and level).
#' @param ... unused.
#' @return a tibble with `response`, `shock` (factors in variable order),
#'   `horizon` and `estimate`; with bands, additionally `level`, `conf.low`,
#'   `conf.high` (one row per cell and level).
#' @family tidiers
#' @examples
#' r <- irf(varm(yraw, p = 2), horizon = 4, ci = "mc", nboot = 50)
#' tidy(r)
#' tidy(r, bands = FALSE)
#' @name tidy.varm_irf
NULL

# Stack percentile bands (last array dimension = sorted probabilities) as rows
# with `level`, `conf.low`, `conf.high`; `est` has one row per array cell.
add_levels <- function(est, bands) {
  nd <- length(dim(bands))
  probs <- as.numeric(dimnames(bands)[[nd]])
  n <- length(probs) / 2
  slice <- function(b) {
    as.vector(do.call(`[`, c(list(bands), rep(list(TRUE), nd - 1L), b)))
  }
  pieces <- lapply(seq_len(n), function(b) {
    est$level <- band_levels(probs)[b]
    est$conf.low <- slice(b)
    est$conf.high <- slice(2 * n + 1 - b)
    est
  })
  do.call(rbind, pieces)
}

tidy_dyn <- function(x, bands = TRUE, ...) {
  vn <- x$varnames
  grid <- expand.grid(
    response = factor(vn, vn),
    shock = factor(vn, vn),
    horizon = 0:x$H
  )
  est <- tibble::as_tibble(cbind(grid, estimate = as.vector(x$estimate)))
  if (!bands || is.null(x$bands)) {
    return(est)
  }
  add_levels(est, x$bands)
}

#' @rdname tidy.varm_irf
#' @export
tidy.varm_irf <- tidy_dyn

#' @rdname tidy.varm_irf
#' @export
tidy.varm_fevd <- tidy_dyn

#' @rdname tidy.varm_irf
#' @param row.names,optional unused.
#' @export
as.data.frame.varm_dyn <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(tidy_dyn(x, ...))
}

#' @exportS3Method tibble::as_tibble
as_tibble.varm_irf <- function(x, ...) tidy_dyn(x, ...)

#' @exportS3Method tibble::as_tibble
as_tibble.varm_fevd <- function(x, ...) tidy_dyn(x, ...)

# historical decomposition -----------------------------------------------------

#' Tidy historical decomposition
#'
#' @param x a `varm_hd` object from [hd()].
#' @param ... unused.
#' @return a long tibble with `time`, `response`, `shock`, `contribution` and
#'   `actual`. The deterministic/initial-condition part appears as
#'   `shock == "baseline"` rows (when available).
#' @family tidiers
#' @examples
#' tidy(hd(varm(yraw, p = 2)))
#' @export
tidy.varm_hd <- function(x, ...) {
  vn <- x$varnames
  Tn <- nrow(x$actual)
  time <- x$dates %||% seq_len(Tn)
  long <- function(shock, value) {
    tibble::tibble(
      time = rep(time, times = length(vn)),
      response = factor(rep(vn, each = Tn), vn),
      shock = shock,
      contribution = as.vector(value),
      actual = as.vector(x$actual)
    )
  }
  parts <- lapply(seq_along(vn), function(j) {
    long(factor(vn[j], c(vn, "baseline")), x$contributions[,, j])
  })
  if (!is.null(x$baseline)) {
    parts <- c(
      parts,
      list(long(factor("baseline", c(vn, "baseline")), x$baseline))
    )
  }
  do.call(rbind, parts)
}

# forecasts -------------------------------------------------------------------

#' Tidy forecasts
#'
#' @param x a `varm_forecast` object from [forecast()].
#' @param ... unused.
#' @return a tibble with `horizon`, `variable`, `estimate` (the mean forecast)
#'   and, when intervals exist, `level`, `conf.low`, `conf.high` (one row per
#'   cell and level).
#' @family tidiers
#' @examples
#' tidy(forecast(varm(yraw, p = 2), h = 4, nboot = 100))
#' @export
tidy.varm_forecast <- function(x, ...) {
  vn <- x$varnames
  h <- nrow(x$mean)
  est <- tibble::tibble(
    horizon = rep(seq_len(h), times = length(vn)),
    variable = factor(rep(vn, each = h), vn),
    estimate = as.vector(x$mean)
  )
  if (is.null(x$bands)) {
    return(est)
  }
  add_levels(est, x$bands)
}

#' @exportS3Method tibble::as_tibble
as_tibble.varm_forecast <- function(x, ...) tidy.varm_forecast(x, ...)

#' @exportS3Method tibble::as_tibble
as_tibble.varm_hd <- function(x, ...) tidy.varm_hd(x, ...)
