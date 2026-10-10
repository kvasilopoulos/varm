#' Tidy, glance at and plot a Bayesian VAR
#'
#' `tidy()` summarises the posterior of the coefficients: `response`, `term`,
#' `estimate` (posterior mean), `std.error` (posterior sd), `conf.low`,
#' `conf.high` (equal-tailed interval). `glance()` is a one-row summary:
#' `K`, `p`, `nobs`, `ndraw`, `method`, `logml`, `accept`, plus the posterior
#' means of the hyperparameters when they were estimated.
#'
#' `autoplot()` shows the coefficient intervals (`type = "coef"`), the posterior
#' densities of the hyperparameters (`"hyper"`) or their trace plots
#' (`"trace"`; the coefficient draws if there are no hyperparameters).
#'
#' @param x,object a [bvarm()] fit.
#' @param conf.level credible level of the interval (default 0.90).
#' @param type `"coef"`, `"hyper"` or `"trace"`.
#' @param ... unused (`plot()` passes them to `autoplot()`).
#' @return `tidy()`, `glance()`: a tibble; `autoplot()`: a `ggplot`.
#' @family Bayesian VAR
#' @importFrom rlang .data
#' @examples
#' \donttest{
#' fit <- bvarm(yraw, p = 2, ndraw = 500, nburn = 300)
#' tidy(fit)
#' glance(fit)
#' autoplot(fit, "hyper")
#' }
#' @name tidy.bvarm
NULL

#' @rdname tidy.bvarm
#' @exportS3Method generics::tidy
tidy.bvarm <- function(x, conf.level = 0.90, ...) {
  B <- x$draws$B
  rn <- rownames(B)
  cn <- colnames(B)
  flat <- matrix(B, dim(B)[1] * dim(B)[2])
  q <- apply(
    flat,
    1,
    stats::quantile,
    probs = c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  )
  tibble::tibble(
    response = rep(cn, each = length(rn)),
    term = rep(rn, times = length(cn)),
    estimate = rowMeans(flat),
    std.error = apply(flat, 1, stats::sd),
    conf.low = q[1, ],
    conf.high = q[2, ]
  )
}

#' @rdname tidy.bvarm
#' @exportS3Method tibble::as_tibble
as_tibble.bvarm <- function(x, ...) tidy.bvarm(x, ...)

#' @rdname tidy.bvarm
#' @exportS3Method generics::glance
glance.bvarm <- function(x, ...) {
  out <- tibble::tibble(
    K = x$K,
    p = x$p,
    nobs = x$nobs,
    ndraw = dim(x$draws$B)[3],
    method = x$draws$method,
    logml = if (is.null(x$logml)) NA_real_ else x$logml,
    accept = if (is.null(x$accept)) NA_real_ else as.numeric(x$accept)
  )
  if (is.null(x$hyper)) {
    return(out)
  }
  tibble::as_tibble(c(out, as.list(colMeans(x$hyper))))
}

#' @rdname tidy.bvarm
#' @exportS3Method ggplot2::autoplot
autoplot.bvarm <- function(
  object,
  type = c("coef", "hyper", "trace"),
  conf.level = 0.90,
  ...
) {
  type <- match.arg(type)
  if (type == "coef") {
    return(bvarm_coef_plot(tidy.bvarm(object, conf.level)))
  }
  if (type == "hyper" && is.null(object$hyper)) {
    rlang::abort("no hyperparameter draws: the prior is fixed or triangular")
  }
  a <- as_chain_array(object)
  d <- chain_long(a[,, seq_len(min(12L, dim(a)[3])), drop = FALSE])
  if (type == "trace") {
    return(trace_plot(d))
  }
  ggplot2::ggplot(d, ggplot2::aes(.data$value)) +
    ggplot2::geom_density(fill = varm_palette(1), alpha = 0.3, colour = NA) +
    ggplot2::facet_wrap("parameter", scales = "free") +
    ggplot2::labs(x = NULL, y = "Posterior density") +
    theme_varm()
}

bvarm_coef_plot <- function(d) {
  d$term <- factor(d$term, levels = rev(unique(d$term)))
  d$response <- factor(d$response, levels = unique(d$response))
  ggplot2::ggplot(d, ggplot2::aes(.data$estimate, .data$term)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high),
      colour = varm_palette(1),
      size = 0.3
    ) +
    ggplot2::facet_wrap("response") +
    ggplot2::labs(x = "Posterior mean and credible interval", y = NULL) +
    theme_varm()
}

#' @rdname tidy.bvarm
#' @export
plot.bvarm <- function(x, ...) plot_gg(x, ...)
