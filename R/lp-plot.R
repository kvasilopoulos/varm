# ggplot2 graphics for local-projection results. Every autoplot() consumes tidy().

#' Plot local-projection results
#'
#' `autoplot()` returns a `ggplot`; `plot()` prints it and returns it invisibly.
#' Shaded ribbons are the stored band levels (widest lightest).
#'
#' * `lpm`, `lp_joint`, `lp_panel`: response paths, one panel per response.
#' * `lp_state`: the high and low state overlaid (`type = "states"`) or their
#'   difference (`type = "diff"`).
#' * `lp_smooth`: smoothed against raw local projection.
#' * `lp_supt`: pointwise against sup-t simultaneous bands.
#' * `lp_forecast`: direct forecast fan chart.
#'
#' @param object,x a result of [lp()], [lp_state()], [lp_smooth()], [lp_panel()],
#'   [lp_joint()], [lp_supt()] or [lp_forecast()].
#' @param type `lp_state` only: `"states"` or `"diff"`.
#' @param ... unused (`plot()` passes them to `autoplot()`).
#' @return a `ggplot` object.
#' @family local projections
#' @importFrom rlang .data
#' @examples
#' fit <- lp(yraw, shock = rate, horizon = 8, lags = 2, contemp = c(infl, unemp))
#' autoplot(fit)
#' @name autoplot.lpm
NULL

# Nested ribbons, one geom per band level (widest first); `fill` optionally names
# a column mapped to the fill aesthetic.
layer_bands <- function(data, fill = NULL, alpha = 0.25) {
  if (!"level" %in% names(data)) {
    data$level <- 0
  }
  data <- data[!is.na(data$level), ]
  one <- function(l) {
    d <- data[data$level == l, ]
    bounds <- ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high)
    if (is.null(fill)) {
      return(ggplot2::geom_ribbon(
        bounds,
        d,
        fill = varm_palette(1),
        alpha = alpha,
        colour = NA
      ))
    }
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$conf.low,
        ymax = .data$conf.high,
        fill = .data[[fill]]
      ),
      d,
      alpha = alpha,
      colour = NA
    )
  }
  lapply(sort(unique(data$level), decreasing = TRUE), one)
}

# Horizon on x, estimate on y, panels by `facet`
lp_paths_plot <- function(data, shock, facet = "response", colour = NULL) {
  aes_line <- if (is.null(colour)) {
    ggplot2::aes(.data$horizon, .data$estimate)
  } else {
    ggplot2::aes(.data$horizon, .data$estimate, colour = .data[[colour]])
  }
  p <- ggplot2::ggplot(data, aes_line) +
    layer_bands(data, fill = colour) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::facet_wrap(facet, scales = "free_y") +
    ggplot2::labs(x = "Horizon (periods)", y = paste("Response to", shock)) +
    theme_varm()
  if (is.null(colour)) {
    return(p + ggplot2::geom_line(colour = varm_palette(1), linewidth = 0.8))
  }
  p + ggplot2::geom_line(linewidth = 0.8)
}

#' @rdname autoplot.lpm
#' @exportS3Method ggplot2::autoplot
autoplot.lpm <- function(object, ...) {
  lp_paths_plot(tidy.lpm(object), object$shock)
}

#' @rdname autoplot.lpm
#' @exportS3Method ggplot2::autoplot
autoplot.lp_state <- function(object, type = c("states", "diff"), ...) {
  type <- match.arg(type)
  d <- tidy.lp_state(object)
  pal <- stats::setNames(varm_palette(2), c("high", "low"))
  if (type == "diff") {
    d <- d[d$state == "diff", ]
    return(
      lp_paths_plot(d, object$shock) +
        ggplot2::labs(
          y = paste("Difference low - high, response to", object$shock)
        )
    )
  }
  d <- d[d$state != "diff", ]
  d$state <- droplevels(d$state)
  lp_paths_plot(d, object$shock, colour = "state") +
    ggplot2::scale_colour_manual(values = pal, name = "State") +
    ggplot2::scale_fill_manual(values = pal, name = "State")
}

#' @rdname autoplot.lpm
#' @exportS3Method ggplot2::autoplot
autoplot.lp_smooth <- function(object, ...) {
  d <- tidy.lp_smooth(object, raw = TRUE)
  pal <- c(smooth = varm_palette(1), raw = "grey40")
  ggplot2::ggplot(d, ggplot2::aes(.data$horizon, .data$estimate)) +
    layer_bands(d[d$type == "smooth", ]) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_line(
      ggplot2::aes(colour = .data$type, linetype = .data$type),
      linewidth = 0.8
    ) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c(smooth = "solid", raw = "dashed"),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Horizon (periods)",
      y = paste("Response to", object$shock)
    ) +
    theme_varm()
}

#' @rdname autoplot.lpm
#' @exportS3Method ggplot2::autoplot
autoplot.lp_supt <- function(object, ...) {
  d <- tidy.lp_supt(object)
  pal <- stats::setNames(varm_palette(2), c("sup-t", "pointwise"))
  ggplot2::ggplot(d, ggplot2::aes(.data$horizon, .data$estimate)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$conf.low,
        ymax = .data$conf.high,
        fill = .data$band
      ),
      alpha = 0.3
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap("response", scales = "free_y") +
    ggplot2::scale_fill_manual(values = pal, name = "Band") +
    ggplot2::labs(x = "Horizon (periods)", y = "Response") +
    theme_varm()
}

#' @rdname autoplot.lpm
#' @exportS3Method ggplot2::autoplot
autoplot.lp_forecast <- function(object, ...) {
  d <- tidy.lp_forecast(object)
  ggplot2::ggplot(d, ggplot2::aes(.data$horizon, .data$estimate)) +
    layer_bands(d) +
    ggplot2::geom_line(colour = varm_palette(1), linewidth = 0.8) +
    ggplot2::facet_wrap("variable", scales = "free_y") +
    ggplot2::labs(x = "Horizon (periods)", y = "Direct forecast") +
    theme_varm()
}

#' @rdname autoplot.lpm
#' @export
plot.lpm <- function(x, ...) plot_gg(x, ...)

#' @rdname autoplot.lpm
#' @export
plot.lp_state <- function(x, ...) plot_gg(x, ...)

#' @rdname autoplot.lpm
#' @export
plot.lp_supt <- function(x, ...) plot_gg(x, ...)

#' @rdname autoplot.lpm
#' @export
plot.lp_forecast <- function(x, ...) plot_gg(x, ...)
