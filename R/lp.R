#' Local projections
#'
#' Jorda (2005) local projections: for each horizon \eqn{h=0,\dots,H} and each response
#' \eqn{y_k},
#' \deqn{y_{k,t+h} = \alpha_h + \beta_h s_t + \gamma_h' w_t + e_{t+h},}
#' where \eqn{s_t} is the shock, \eqn{w_t} the controls (lags of `Y`, optional contemporaneous
#' variables, extra controls). \eqn{\beta_h} is the impulse response at horizon \eqn{h}.
#' The regressions, 2SLS and HAC standard errors run in C++.
#'
#' @param Y numeric matrix / data frame (T x K) of responses (and the source of lag controls).
#' @param shock the shock \eqn{s_t}: a column name/index of `Y` (observed shock; unit-impact
#'   normalisation) or a numeric vector of length T (external shock series, e.g. narrative).
#' @param horizon maximum horizon H.
#' @param lags number of lags of `Y` in the controls (and of an external shock).
#' @param contemp names/indices of `Y` columns entered contemporaneously as controls. For a
#'   recursive (Cholesky) identification with `shock` a column of `Y`, list the variables ordered
#'   *before* the shock.
#' @param controls optional extra T x q matrix of controls.
#' @param iv optional T x q matrix of instruments for `shock` (LP-IV; Stock & Watson 2018).
#' @param responses subset of `Y` columns to use as responses (default all).
#' @param cumulative variables whose response is cumulated (\eqn{\sum_{j=0}^h y_{t+j}}), for
#'   variables measured in differences.
#' @param se `"nw"` Newey-West/Bartlett or `"hc"` Eicker-Huber-White (HC1).
#' @param nw_lag Newey-West truncation lag; default `h + 1`.
#' @param lag_augment Montiel Olea & Plagborg-Moller (2021) lag augmentation: add one extra lag
#'   that is not interpreted and use `"hc"` standard errors, valid without stationarity
#'   assumptions and robust at long horizons.
#' @param level confidence levels for the pointwise bands.
#' @param dof_adjust apply an `n/(n-m)` small-sample factor to the Newey-West variance.
#' @return an `lpm` object: `$estimate` (`(H+1) x K`), `$se`, `$bands`, `$nobs`.
#' @examples
#' fit <- lp(yraw, shock = "V3", horizon = 12, lags = 2, contemp = c("V1", "V2"))
#' fit
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
  cl <- match.call()
  se <- match.arg(se)
  Y <- as_data_matrix(Y)
  n <- nrow(Y)
  vn <- colnames(Y)
  if (lag_augment) {
    lags <- lags + 1L
    se <- "hc"
  }
  # shock
  if (length(shock) == 1L && (is.character(shock) || shock %% 1 == 0)) {
    sj <- resolve_var(shock, vn, "shock")
    s <- Y[, sj]
    sname <- vn[sj]
    external <- FALSE
  } else {
    s <- as.numeric(shock)
    if (length(s) != n) {
      stop("external `shock` must have length nrow(Y)", call. = FALSE)
    }
    sname <- "shock"
    external <- TRUE
  }
  if (anyNA(s)) {
    stop("`shock` contains missing values", call. = FALSE)
  }
  lagY <- do.call(
    cbind,
    lapply(seq_len(lags), function(l) {
      rbind(matrix(NA, l, ncol(Y)), Y[seq_len(n - l), , drop = FALSE])
    })
  )
  W <- cbind(const = 1, lagY)
  colnames(W) <- c(
    "const",
    unlist(lapply(seq_len(lags), function(l) paste0(vn, ".l", l)))
  )
  if (external) {
    ls_ <- do.call(
      cbind,
      lapply(seq_len(lags), function(l) c(rep(NA, l), s[seq_len(n - l)]))
    )
    colnames(ls_) <- paste0("shock.l", seq_len(lags))
    W <- cbind(W, ls_)
  }
  if (!is.null(contemp)) {
    cj <- resolve_var(contemp, vn, "contemporaneous variable")
    W <- cbind(W, Y[, cj, drop = FALSE])
  }
  if (!is.null(controls)) {
    W <- cbind(W, as_data_matrix(controls, "controls"))
  }
  Z <- if (is.null(iv)) matrix(0, n, 0L) else as_data_matrix(iv, "iv")
  rs <- if (is.null(responses)) {
    seq_len(ncol(Y))
  } else {
    resolve_var(responses, vn, "response")
  }
  cum <- as.integer(
    rs %in%
      if (is.null(cumulative)) {
        integer(0)
      } else {
        resolve_var(cumulative, vn, "cumulative variable")
      }
  )
  H <- as.integer(horizon)
  r <- cpp_lp(
    Y[, rs, drop = FALSE],
    matrix(s, ncol = 1L),
    W,
    Z,
    H,
    as.integer(lags),
    cum,
    if (se == "nw") 0L else 1L,
    if (is.null(nw_lag)) -1L else as.integer(nw_lag),
    dof_adjust
  )
  est <- r$beta[,, 1]
  sev <- r$se[,, 1]
  est <- matrix(est, H + 1L, length(rs), dimnames = list(0:H, vn[rs]))
  sev <- matrix(sev, H + 1L, length(rs), dimnames = dimnames(est))
  probs <- level_probs(level)
  z <- stats::qnorm(probs)
  bands <- array(
    NA_real_,
    c(H + 1L, length(rs), length(probs)),
    dimnames = list(0:H, vn[rs], format(probs))
  )
  for (i in seq_along(z)) {
    bands[,, i] <- est + z[i] * sev
  }
  structure(
    list(
      estimate = est,
      se = sev,
      bands = bands,
      probs = probs,
      level = level,
      nobs = r$nobs,
      H = H,
      shock = sname,
      responses = vn[rs],
      lags = lags,
      se_type = se,
      iv = !is.null(iv),
      call = cl
    ),
    class = "lpm"
  )
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
  print(round(
    x$estimate[unique(pmin(x$H + 1L, c(1, 2, 5, 9, 13, 21))), , drop = FALSE],
    digits
  ))
  invisible(x)
}

#' @export
plot.lpm <- function(x, ...) {
  k <- length(x$responses)
  hh <- 0:x$H
  op <- graphics::par(
    mfrow = c(1, k),
    mar = c(2.4, 2.6, 1.6, 0.6),
    mgp = c(1.4, 0.4, 0)
  )
  on.exit(graphics::par(op))
  nb <- dim(x$bands)[3] / 2
  for (j in seq_len(k)) {
    graphics::plot(
      hh,
      x$estimate[, j],
      type = "n",
      ylim = range(x$bands[, j, ], x$estimate[, j]),
      main = paste0(x$shock, " -> ", x$responses[j]),
      xlab = "",
      ylab = "",
      cex.main = 0.85
    )
    for (b in seq_len(nb)) {
      graphics::polygon(
        c(hh, rev(hh)),
        c(x$bands[, j, b], rev(x$bands[, j, dim(x$bands)[3] + 1 - b])),
        border = NA,
        col = grDevices::adjustcolor("firebrick", 0.22)
      )
    }
    graphics::abline(h = 0, col = "grey50")
    graphics::lines(hh, x$estimate[, j], lwd = 2, col = "firebrick4")
  }
  invisible(x)
}

#' @export
as.data.frame.lpm <- function(x, row.names = NULL, optional = FALSE, ...) {
  g <- expand.grid(h = 0:x$H, response = x$responses, stringsAsFactors = FALSE)
  df <- data.frame(g, estimate = as.vector(x$estimate), se = as.vector(x$se))
  for (b in seq_along(x$probs)) {
    df[[paste0("q", format(x$probs[b]))]] <- as.vector(x$bands[,, b])
  }
  df
}
