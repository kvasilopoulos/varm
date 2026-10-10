#' Tidy and glance at local-projection results
#'
#' `tidy()` returns one row per horizon, response and (for [lp_state()]) state,
#' repeated for every confidence level so that bands are tidy too:
#' `horizon`, `response`, (`state`), `estimate`, `std.error`, `level`,
#' `conf.low`, `conf.high`. `glance()` returns a one-row summary of the fit.
#' `as_tibble()` is an alias of `tidy()`.
#'
#' @param x a result of [lp()], [lp_joint()], [lp_state()], [lp_panel()],
#'   [lp_smooth()], [lp_supt()] or [lp_forecast()].
#' @param conf.int include the `level`, `conf.low` and `conf.high` columns
#'   (one row per band level).
#' @param raw `lp_smooth` only: also return the unsmoothed local projection as
#'   rows with `type == "raw"` (no bands).
#' @param ... unused.
#' @return a tibble.
#' @family local projections
#' @examples
#' fit <- lp(yraw, shock = rate, horizon = 6, lags = 2, contemp = c(infl, unemp))
#' tidy(fit)
#' glance(fit)
#' @name tidy.lpm
NULL

# Long tibble of an estimate array (H+1 x K [x S]) with its band array
tidy_cells <- function(est, se, bands = NULL, probs = NULL) {
  dn <- dimnames(est)
  H1 <- nrow(est)
  K <- ncol(est)
  S <- length(est) / (H1 * K)
  states <- if (length(dim(est)) == 3L) dn[[3]]
  out <- tibble::tibble(
    horizon = rep(as.integer(dn[[1]]), K * S),
    response = factor(rep(rep(dn[[2]], each = H1), S), levels = dn[[2]]),
    state = if (!is.null(states)) {
      factor(rep(states, each = H1 * K), levels = states)
    },
    estimate = as.vector(est),
    std.error = as.vector(se)
  )
  if (is.null(bands)) {
    return(out)
  }
  np <- length(probs)
  bl <- matrix(bands, ncol = np)
  one <- function(b) {
    o <- out
    o$level <- round(1 - 2 * probs[b], 8)
    o$conf.low <- bl[, b]
    o$conf.high <- bl[, np + 1L - b]
    o
  }
  tibble::as_tibble(do.call(rbind, lapply(seq_len(np %/% 2L), one)))
}

#' @rdname tidy.lpm
#' @exportS3Method generics::tidy
tidy.lpm <- function(x, conf.int = TRUE, ...) {
  tidy_cells(x$estimate, x$se, if (conf.int) x$bands, x$probs)
}

#' @rdname tidy.lpm
#' @exportS3Method generics::tidy
tidy.lp_state <- function(x, conf.int = TRUE, ...) {
  tidy_cells(x$estimate, x$se, if (conf.int) x$bands, x$probs)
}

#' @rdname tidy.lpm
#' @exportS3Method generics::tidy
tidy.lp_smooth <- function(x, conf.int = TRUE, raw = FALSE, ...) {
  out <- tidy_cells(x$estimate, x$se, if (conf.int) x$bands, x$probs)
  if (!raw) {
    return(out)
  }
  out$type <- "smooth"
  unsmoothed <- tibble::tibble(
    horizon = 0:x$H,
    response = out$response[1],
    estimate = x$lp,
    std.error = x$lp_se,
    level = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    type = "raw"
  )
  tibble::as_tibble(rbind(out, unsmoothed[names(out)]))
}

#' @rdname tidy.lpm
#' @exportS3Method generics::tidy
tidy.lp_forecast <- function(x, conf.int = TRUE, ...) {
  out <- tidy_cells(x$mean, x$se, if (conf.int) x$bands, x$probs)
  names(out)[names(out) == "response"] <- "variable"
  out
}

#' @rdname tidy.lpm
#' @exportS3Method generics::tidy
tidy.lp_supt <- function(x, ...) {
  df <- tibble::as_tibble(as.data.frame(x))
  key <- c("response", "horizon", "estimate", "std.error")
  band <- function(name, lo, hi) {
    out <- df[key]
    out$band <- name
    out$level <- attr(x, "level")
    out$conf.low <- df[[lo]]
    out$conf.high <- df[[hi]]
    out
  }
  out <- rbind(
    band("pointwise", "pw_low", "pw_high"),
    band("sup-t", "conf.low", "conf.high")
  )
  out$band <- factor(out$band, c("sup-t", "pointwise"))
  tibble::as_tibble(out)
}

#' @rdname tidy.lpm
#' @exportS3Method tibble::as_tibble
as_tibble.lpm <- function(x, ...) tidy.lpm(x, ...)

#' @rdname tidy.lpm
#' @exportS3Method tibble::as_tibble
as_tibble.lp_state <- function(x, ...) tidy.lp_state(x, ...)

#' @rdname tidy.lpm
#' @exportS3Method tibble::as_tibble
as_tibble.lp_forecast <- function(x, ...) tidy.lp_forecast(x, ...)

#' @export
as.data.frame.lpm <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(tidy.lpm(x, ...))
}

#' @rdname tidy.lpm
#' @exportS3Method generics::glance
glance.lpm <- function(x, ...) {
  tibble::tibble(
    shock = x$shock,
    horizon = x$H,
    lags = x$lags,
    se_type = x$se_type,
    iv = x$iv,
    nobs = max(x$nobs)
  )
}

#' @rdname tidy.lpm
#' @exportS3Method generics::glance
glance.lp_panel <- function(x, ...) {
  tibble::tibble(
    shock = x$shock,
    horizon = x$H,
    fe = x$fe,
    cluster = x$cluster,
    nobs = x$nobs[1],
    nclusters = x$nclusters[1]
  )
}

#' @rdname tidy.lpm
#' @exportS3Method generics::glance
glance.lp_smooth <- function(x, ...) {
  tibble::tibble(
    shock = x$shock,
    horizon = x$H,
    nobs = x$nobs[1],
    lambda = x$lambda,
    edf = x$edf
  )
}

#' @rdname tidy.lpm
#' @exportS3Method generics::glance
glance.lp_state <- function(x, ...) {
  tibble::tibble(
    shock = x$shock,
    horizon = x$H,
    lags = x$lags,
    timing = x$timing,
    transition = x$transition,
    cutoff = x$cutoff,
    share_low = x$share_low,
    nobs = max(x$nobs),
    wald_p = min(x$wald[, "p"])
  )
}
