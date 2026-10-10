#' Tidy and plot a TVP-VAR
#'
#' `tidy()` summarises the posterior paths through time as a tibble:
#' `type = "coef"` gives `date`, `response` (equation), `term`, `estimate`
#' (posterior median), `conf.low`, `conf.high`; `"vol"` and `"shock_vol"` give
#' `date`, `variable`, `estimate`, `conf.low`, `conf.high` for the reduced-form
#' residual and the orthogonalised shock standard deviations
#' ([tvp_volatility()]). `glance()` is a one-row summary.
#'
#' `autoplot()` draws the paths with their credible band (`type = "coef"`,
#' `"vol"`, `"shock_vol"`) or the response of `response` to `shock` through time
#' at the given `horizon`s (`type = "irf"`, via [irf_path()]).
#'
#' @param x,object a [tvpvar()] fit.
#' @param type `"coef"`, `"vol"`, `"shock_vol"` (and `"irf"` for `autoplot()`).
#' @param equation <[`tidy-select`][tidyselect::language]> equation(s) for
#'   `type = "coef"`. `tidy()` defaults to all, `autoplot()` to the first.
#' @param level credible level of the band.
#' @param const include the intercept path (`type = "coef"`).
#' @param shock,response <[`tidy-select`][tidyselect::language]> variables of the
#'   response shown by `type = "irf"`.
#' @param horizon horizons of the response shown by `type = "irf"`.
#' @param ... unused (`plot()` passes them to `autoplot()`).
#' @return `tidy()`, `glance()`: a tibble; `autoplot()`: a `ggplot`.
#' @family TVP-VAR
#' @importFrom rlang .data
#' @examples
#' \donttest{
#' fit <- tvpvar(yraw, p = 2, ndraw = 200, nburn = 200)
#' tidy(fit, "vol")
#' autoplot(fit, "coef", equation = infl)
#' }
#' @name tidy.tvpvar
NULL

#' @rdname tidy.tvpvar
#' @exportS3Method generics::tidy
tidy.tvpvar <- function(
  x,
  type = c("coef", "vol", "shock_vol"),
  equation = NULL,
  level = 0.68,
  const = TRUE,
  ...
) {
  type <- match.arg(type)
  if (type != "coef") {
    return(tvp_volatility(
      x,
      if (type == "vol") "reduced" else "structural",
      level
    ))
  }
  eq <- select_or(rlang::enquo(equation), x$varnames, "equation", seq_len(x$K))
  terms <- if (const || !x$const) seq_len(x$m) else seq_len(x$m - 1L)
  q <- tvp_quantiles(x$beta[terms, eq, , , drop = FALSE], level)
  nt <- length(terms)
  ne <- length(eq)
  tibble::tibble(
    date = rep(x$dates, each = nt * ne),
    response = factor(
      rep(rep(x$varnames[eq], each = nt), x$T),
      levels = x$varnames[eq]
    ),
    term = factor(
      rep(x$coefnames[terms], ne * x$T),
      levels = x$coefnames[terms]
    ),
    estimate = as.vector(q[2, , , ]),
    conf.low = as.vector(q[1, , , ]),
    conf.high = as.vector(q[3, , , ])
  )
}

#' @rdname tidy.tvpvar
#' @exportS3Method tibble::as_tibble
as_tibble.tvpvar <- function(x, ...) tidy.tvpvar(x, ...)

#' @rdname tidy.tvpvar
#' @exportS3Method generics::glance
glance.tvpvar <- function(x, ...) {
  tibble::tibble(
    K = x$K,
    p = x$p,
    nobs = x$T,
    ndraw = x$ndraw,
    npaths = dim(x$h)[3],
    tau = x$tau,
    tv_beta = x$tv_beta,
    tv_A = x$tv_A,
    sv = x$sv
  )
}

# Posterior median and band against date, one panel per `facet`
tvp_paths_plot <- function(d, facet, ylab, colour = NULL) {
  aes_line <- ggplot2::aes(.data$date, .data$estimate)
  if (!is.null(colour)) {
    aes_line <- ggplot2::aes(
      .data$date,
      .data$estimate,
      colour = .data[[colour]]
    )
  }
  p <- ggplot2::ggplot(d, aes_line) +
    layer_bands(d, fill = colour) +
    ggplot2::labs(x = NULL, y = ylab) +
    theme_varm()
  p <- p +
    if (is.null(colour)) {
      ggplot2::geom_line(colour = varm_palette(1))
    } else {
      ggplot2::geom_line()
    }
  if (is.null(facet)) {
    return(p)
  }
  p + facet
}

#' @rdname tidy.tvpvar
#' @exportS3Method ggplot2::autoplot
autoplot.tvpvar <- function(
  object,
  type = c("coef", "vol", "shock_vol", "irf"),
  equation = NULL,
  level = 0.68,
  const = FALSE,
  shock = NULL,
  response = NULL,
  horizon = c(0L, 4L, 8L),
  ...
) {
  type <- match.arg(type)
  if (type == "irf") {
    return(tvp_irf_plot(
      object,
      rlang::enquo(shock),
      rlang::enquo(response),
      horizon,
      level
    ))
  }
  if (type != "coef") {
    d <- tidy.tvpvar(object, type, level = level)
    return(tvp_paths_plot(
      d,
      ggplot2::facet_wrap("variable", scales = "free_y"),
      "Standard deviation"
    ))
  }
  eq <- rlang::enquo(equation)
  if (rlang::quo_is_null(eq)) {
    eq <- rlang::quo(1L)
  }
  d <- tidy.tvpvar(object, "coef", !!eq, level, const)
  facet <- if (nlevels(d$response) > 1L) {
    ggplot2::facet_grid(response ~ term, scales = "free_y")
  } else {
    ggplot2::facet_wrap("term", scales = "free_y")
  }
  tvp_paths_plot(
    d,
    facet,
    paste("Equation:", paste(levels(d$response), collapse = ", "))
  )
}

tvp_irf_plot <- function(object, shock, response, horizon, level) {
  if (rlang::quo_is_null(shock) || rlang::quo_is_null(response)) {
    rlang::abort("`type = \"irf\"` needs `shock` and `response`")
  }
  paths <- lapply(horizon, function(h) {
    irf_path(object, !!shock, !!response, horizon = h, level = level)
  })
  d <- do.call(
    rbind,
    lapply(seq_along(horizon), function(i) {
      out <- tibble::as_tibble(paths[[i]])
      out$horizon <- factor(horizon[i], levels = horizon)
      out
    })
  )
  label <- sub(", h = .*", "", attr(paths[[1]], "label"))
  pal <- stats::setNames(varm_palette(length(horizon)), horizon)
  p <- tvp_paths_plot(d, NULL, label, colour = "horizon")
  p +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::scale_colour_manual(values = pal, name = "Horizon") +
    ggplot2::scale_fill_manual(values = pal, name = "Horizon")
}

#' @rdname tidy.tvpvar
#' @exportS3Method ggplot2::autoplot
autoplot.tvpvar_path <- function(object, ...) {
  d <- tibble::as_tibble(object)
  tvp_paths_plot(d, NULL, attr(object, "label")) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50")
}

#' @rdname tidy.tvpvar
#' @export
plot.tvpvar <- function(x, ...) plot_gg(x, ...)

#' @rdname tidy.tvpvar
#' @export
plot.tvpvar_path <- function(x, ...) plot_gg(x, ...)
