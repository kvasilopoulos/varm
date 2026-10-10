# ggplot2 graphics: every plot is built from the tidy() output -----------------

#' Plot impulse responses and variance decompositions
#'
#' Grid of panels (rows: responding variable, columns: shock) with the point or
#' median response and nested ribbons for each stored level, widest lightest.
#' `plot()` prints the same ggplot and returns it invisibly.
#'
#' @param object,x a `varm_irf` or `varm_fevd` object.
#' @param shocks,responses <[`tidy-select`][tidyselect::language]> subsets of
#'   the shocks and responding variables, e.g. `shocks = rate`.
#' @param ci `"band"` (shaded ribbons, default), `"line"` (dashed bounds) or
#'   `"none"`.
#' @param compare named list of further `varm_irf` / `varm_fevd` objects to
#'   overlay on `object` (labelled "Baseline" and the list names), one colour
#'   each.
#' @param ... unused.
#' @return A `ggplot` object.
#' @family plotting
#' @examples
#' fit <- varm(yraw, p = 2)
#' r <- irf(fit, horizon = 8, ci = "mc", nboot = 100)
#' ggplot2::autoplot(r)
#' ggplot2::autoplot(r, shocks = rate, responses = c(infl, unemp))
#' ggplot2::autoplot(r, ci = "line")
#' r2 <- irf(varm(yraw, p = 4), horizon = 8, ci = "none")
#' ggplot2::autoplot(r, compare = list("p = 4" = r2), ci = "none")
#' ggplot2::autoplot(fevd(fit, horizon = 8, ci = "none"))
#' @name autoplot.varm_irf
NULL

autoplot_dyn <- function(
  object,
  shocks,
  responses,
  ylab,
  zero_line,
  ci = "band",
  compare = NULL
) {
  ci <- match.arg(ci, c("band", "line", "none"))
  vn <- object$varnames
  sh <- vn[select_vars(shocks, vn, "shock") %||% seq_along(vn)]
  rs <- vn[select_vars(responses, vn, "response") %||% seq_along(vn)]
  objs <- c(list(object), compare)
  nms <- names(compare)
  if (length(compare) && (is.null(nms) || anyNA(nms) || any(nms == ""))) {
    rlang::abort("`compare` must be a named list of objects.")
  }
  nms <- c("Baseline", nms)
  d <- do.call(
    rbind,
    Map(
      function(o, n) {
        t <- tidy(o)
        t$model <- n
        t
      },
      objs,
      nms
    )
  )
  d$model <- factor(d$model, nms)
  d <- d[d$shock %in% sh & d$response %in% rs, ]
  d$shock <- factor(d$shock, sh)
  d$response <- factor(d$response, rs)
  multi <- length(objs) > 1L
  has_ci <- ci != "none" && "level" %in% names(d)
  cell <- c("model", "response", "shock", "horizon")
  line <- d[!duplicated(d[cell]), ]
  p <- ggplot2::ggplot(
    line,
    ggplot2::aes(.data$horizon, .data$estimate, colour = .data$model)
  )
  if (zero_line) {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey60")
  }
  if (has_ci && ci == "band") {
    p <- p + band_ribbon(d)
  } else if (has_ci) {
    p <- p + band_lines(d)
  }
  pal <- varm_palette(length(objs))
  p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::guides(colour = if (multi) "legend" else "none") +
    ggplot2::facet_grid(response ~ shock, scales = "free_y") +
    ggplot2::labs(x = "Horizon (periods)", y = ylab) +
    theme_varm()
}

# nested ribbons: one per level (and model), overlapping with transparency
band_ribbon <- function(d) {
  ggplot2::geom_ribbon(
    data = d,
    ggplot2::aes(
      ymin = .data$conf.low,
      ymax = .data$conf.high,
      group = interaction(.data$model, .data$level),
      fill = .data$model
    ),
    alpha = 0.22,
    colour = NA,
    show.legend = FALSE
  )
}

# dashed bound lines instead of shading
band_lines <- function(d) {
  b <- function(y) {
    ggplot2::geom_line(
      data = d,
      ggplot2::aes(
        y = .data[[y]],
        group = interaction(.data$model, .data$level)
      ),
      linetype = "dashed",
      linewidth = 0.4
    )
  }
  list(b("conf.low"), b("conf.high"))
}

#' @rdname autoplot.varm_irf
#' @exportS3Method ggplot2::autoplot
autoplot.varm_irf <- function(
  object,
  shocks = NULL,
  responses = NULL,
  ci = c("band", "line", "none"),
  compare = NULL,
  ...
) {
  autoplot_dyn(
    object,
    rlang::enquo(shocks),
    rlang::enquo(responses),
    "Response",
    TRUE,
    ci,
    compare
  )
}

#' @rdname autoplot.varm_irf
#' @exportS3Method ggplot2::autoplot
autoplot.varm_fevd <- function(
  object,
  shocks = NULL,
  responses = NULL,
  ci = c("band", "line", "none"),
  compare = NULL,
  ...
) {
  autoplot_dyn(
    object,
    rlang::enquo(shocks),
    rlang::enquo(responses),
    "Share of forecast error variance",
    FALSE,
    ci,
    compare
  )
}

# print the autoplot and return it invisibly
plot_gg <- function(x, ...) {
  p <- ggplot2::autoplot(x, ...)
  print(p)
  invisible(p)
}

#' @rdname autoplot.varm_irf
#' @export
plot.varm_dyn <- plot_gg

#' Plot forecasts: history and fan chart
#'
#' The last `n_history` observations are drawn as a line followed by the mean
#' forecast and nested interval ribbons.
#'
#' @param object,x a `varm_forecast` object.
#' @param variables <[`tidy-select`][tidyselect::language]> variables to show.
#' @param n_history number of past observations to show.
#' @param ... unused.
#' @return A `ggplot` object.
#' @family plotting
#' @examples
#' fc <- forecast(varm(yraw, p = 2), h = 8, nboot = 200)
#' ggplot2::autoplot(fc)
#' ggplot2::autoplot(fc, variables = rate, n_history = 20)
#' @exportS3Method ggplot2::autoplot
autoplot.varm_forecast <- function(
  object,
  variables = NULL,
  n_history = 40L,
  ...
) {
  vn <- object$varnames
  keep <- vn[
    select_vars(rlang::enquo(variables), vn, "variable") %||% seq_along(vn)
  ]
  d <- tidy(object)
  d <- d[d$variable %in% keep, ]
  d$variable <- factor(d$variable, keep)
  hist <- history_frame(object, keep, n_history)
  time <- forecast_time(object)
  d$time <- time$future[d$horizon]
  p <- ggplot2::ggplot(
    d[!duplicated(d[c("variable", "horizon")]), ],
    ggplot2::aes(.data$time, .data$estimate)
  )
  if (!is.null(hist)) {
    hist$time <- time$past[hist$row]
    p <- p +
      ggplot2::geom_line(
        data = hist,
        ggplot2::aes(.data$time, .data$value),
        colour = "grey30"
      )
  }
  if ("level" %in% names(d)) {
    d$model <- "Baseline"
    p <- p +
      band_ribbon(d) +
      ggplot2::scale_fill_manual(values = varm_palette(1))
  }
  p +
    ggplot2::geom_line(colour = varm_palette(1), linewidth = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y") +
    ggplot2::labs(x = time$label, y = NULL) +
    theme_varm()
}

# last `n` observations of the chosen variables in long form
history_frame <- function(object, keep, n) {
  Y <- object$history
  if (is.null(Y)) {
    return(NULL)
  }
  rows <- seq.int(max(1L, nrow(Y) - n + 1L), nrow(Y))
  tibble::tibble(
    row = rep(rows, times = length(keep)),
    variable = factor(rep(keep, each = length(rows)), keep),
    value = as.vector(Y[rows, keep, drop = FALSE])
  )
}

# x-axis positions of the observations (`past`) and the forecasts (`future`)
forecast_time <- function(object) {
  Tn <- nrow(object$history %||% object$mean)
  h <- nrow(object$mean)
  dates <- object$dates
  if (inherits(dates, "Date") && length(dates) == Tn) {
    return(list(past = dates, future = future_dates(dates, h), label = "Date"))
  }
  list(past = seq_len(Tn), future = Tn + seq_len(h), label = "Time (periods)")
}

future_dates <- function(dates, h) {
  gap <- stats::median(diff(as.numeric(dates)))
  by <- if (gap >= 365) {
    "year"
  } else if (gap >= 89) {
    "3 months"
  } else if (gap >= 28) {
    "month"
  } else if (gap >= 7) {
    "week"
  } else {
    "day"
  }
  seq(dates[length(dates)], by = by, length.out = h + 1L)[-1L]
}

#' @rdname autoplot.varm_forecast
#' @export
plot.varm_forecast <- plot_gg

#' Plot historical decompositions
#'
#' Stacked columns of the shock contributions (and the baseline) to each
#' variable, with the actual series as a line.
#'
#' @param object,x a `varm_hd` object from [hd()].
#' @param responses <[`tidy-select`][tidyselect::language]> variables to show.
#' @param ... unused.
#' @return A `ggplot` object.
#' @family plotting
#' @examples
#' ggplot2::autoplot(hd(varm(yraw, p = 2)), responses = rate)
#' @exportS3Method ggplot2::autoplot
autoplot.varm_hd <- function(object, responses = NULL, ...) {
  vn <- object$varnames
  keep <- vn[
    select_vars(rlang::enquo(responses), vn, "response") %||% seq_along(vn)
  ]
  d <- tidy(object)
  d <- d[d$response %in% keep, ]
  d$response <- factor(d$response, keep)
  width <- min(diff(as.numeric(sort(unique(d$time)))))
  fills <- c(
    stats::setNames(varm_palette(length(vn)), vn),
    baseline = "grey75"
  )
  ggplot2::ggplot(d, ggplot2::aes(.data$time)) +
    ggplot2::geom_col(
      ggplot2::aes(y = .data$contribution, fill = .data$shock),
      width = width
    ) +
    ggplot2::geom_line(
      data = d[!duplicated(d[c("time", "response")]), ],
      ggplot2::aes(y = .data$actual),
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = fills, name = "Shock") +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$response),
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(x = NULL, y = "Contribution") +
    theme_varm()
}

#' @rdname autoplot.varm_hd
#' @export
plot.varm_hd <- plot_gg

#' Plot parameter draws: trace and density
#'
#' @param object,x a `varm_draws` object from [draws()] (or a Bayesian fit's
#'   draws).
#' @param type `"trace"` or `"density"`.
#' @param terms <[`tidy-select`][tidyselect::language]> coefficient (row) names
#'   such as `infl.l1`, shown for every equation. Default: each variable's own
#'   first lag.
#' @param ... unused.
#' @return A `ggplot` object.
#' @family plotting
#' @examples
#' d <- draws(varm(yraw, p = 2), "mc", n = 200)
#' ggplot2::autoplot(d)
#' ggplot2::autoplot(d, type = "density", terms = c(rate.l1, const))
#' @exportS3Method ggplot2::autoplot
autoplot.varm_draws <- function(
  object,
  type = c("trace", "density"),
  terms = NULL,
  ...
) {
  type <- match.arg(type)
  d <- draws_long(object, rlang::enquo(terms))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$value))
  if (type == "density") {
    p <- p +
      ggplot2::geom_density(
        fill = varm_palette(1),
        alpha = 0.3,
        colour = varm_palette(1)
      ) +
      ggplot2::labs(x = "Value", y = "Density")
  } else {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$draw, y = .data$value),
        colour = varm_palette(1),
        linewidth = 0.3
      ) +
      ggplot2::labs(x = "Draw", y = "Value")
  }
  p +
    ggplot2::facet_wrap(ggplot2::vars(.data$parameter), scales = "free") +
    theme_varm()
}

# long tibble of selected coefficients: draw, parameter, value
draws_long <- function(object, terms) {
  B <- object$B
  K <- object$K
  rn <- rownames(B) %||%
    rownames(object$fit$B) %||%
    paste0("b", seq_len(nrow(B)))
  vn <- object$varnames
  sel <- select_vars(terms, rn, "term")
  pairs <- if (is.null(sel)) {
    cbind(row = seq_len(K), col = seq_len(K))
  } else {
    as.matrix(expand.grid(row = sel, col = seq_len(K)))
  }
  N <- dim(B)[3]
  pieces <- lapply(seq_len(nrow(pairs)), function(k) {
    r <- pairs[k, "row"]
    cl <- pairs[k, "col"]
    tibble::tibble(
      draw = seq_len(N),
      parameter = paste0(vn[cl], ": ", rn[r]),
      value = B[r, cl, ]
    )
  })
  do.call(rbind, pieces)
}

#' @rdname autoplot.varm_draws
#' @export
plot.varm_draws <- plot_gg
