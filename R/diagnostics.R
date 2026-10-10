#' MCMC convergence diagnostics
#'
#' Effective sample size (Geyer's initial monotone positive sequence,
#' multi-chain form of Vehtari et al. 2021), split R-hat, and the Geweke (1992)
#' mean-equality z-score (first 10% vs last 50% of each chain, variances
#' inflated by the autocorrelation time) for every parameter.
#'
#' @param draws one of: a numeric matrix `iterations x parameters` (one chain);
#'   a 3-d array `iterations x chains x parameters`; a list of such matrices
#'   (one per chain); a `bvarm` (hyperparameter draws if hierarchical,
#'   otherwise the coefficient draws `B`); or a `varm_draws` object
#'   (coefficients).
#' @return a tibble of class `mcmc_diagnostics` with columns `term`,
#'   `estimate` (posterior mean), `std.error` (posterior sd), `ess` (bulk
#'   effective sample size over all chains), `rhat` (split R-hat; every chain is
#'   cut in two halves, so a single chain works) and `geweke_z` (mean over
#'   chains, `|z| > 2` suggests non-stationarity). The draws are kept for
#'   [autoplot()][autoplot.mcmc_diagnostics()]; [glance()] summarises the table.
#' @family Bayesian VAR
#' @examples
#' set.seed(1)
#' x <- cbind(a = arima.sim(list(ar = 0.8), 2000), b = rnorm(2000))
#' diag <- mcmc_diagnostics(x)
#' diag
#' glance(diag)
#' @export
mcmc_diagnostics <- function(draws) {
  a <- as_chain_array(draws)
  n <- dim(a)[1]
  m <- dim(a)[2]
  P <- dim(a)[3]
  if (n < 8L) {
    rlang::abort("need at least 8 iterations per chain")
  }
  nm <- dimnames(a)[[3]]
  if (is.null(nm)) {
    nm <- paste0("V", seq_len(P))
  }
  dimnames(a) <- list(NULL, NULL, nm)
  one <- function(k) {
    x <- matrix(a[,, k], n, m)
    tibble::tibble(
      term = nm[k],
      estimate = mean(x),
      std.error = stats::sd(x),
      ess = ess_chains(x),
      rhat = rhat_split(x),
      geweke_z = mean(apply(x, 2, geweke_z), na.rm = TRUE)
    )
  }
  out <- do.call(rbind, lapply(seq_len(P), one))
  structure(
    out,
    draws = a,
    class = c("mcmc_diagnostics", class(out))
  )
}

#' @export
print.mcmc_diagnostics <- function(x, digits = 3, ...) {
  print.data.frame(format(as.data.frame(x), digits = digits), row.names = FALSE)
  invisible(x)
}

#' @rdname mcmc_diagnostics
#' @param x an `mcmc_diagnostics` object.
#' @param ... unused.
#' @exportS3Method generics::tidy
tidy.mcmc_diagnostics <- function(x, ...) {
  tibble::as_tibble(as.data.frame(x))
}

#' @rdname mcmc_diagnostics
#' @exportS3Method tibble::as_tibble
as_tibble.mcmc_diagnostics <- function(x, ...) tidy.mcmc_diagnostics(x, ...)

#' @rdname mcmc_diagnostics
#' @exportS3Method generics::glance
glance.mcmc_diagnostics <- function(x, ...) {
  tibble::tibble(
    n_param = nrow(x),
    n_chain = dim(attr(x, "draws"))[2],
    n_iter = dim(attr(x, "draws"))[1],
    min_ess = min(x$ess, na.rm = TRUE),
    max_rhat = max(x$rhat, na.rm = TRUE),
    max_abs_geweke = max(abs(x$geweke_z), na.rm = TRUE)
  )
}

#' Plot MCMC diagnostics
#'
#' `type = "summary"` shows the effective sample size and R-hat of every
#' parameter, `"trace"` the chains and `"density"` the marginal posterior
#' densities of the selected parameters.
#'
#' @param object an [mcmc_diagnostics()] result.
#' @param type `"summary"`, `"trace"` or `"density"`.
#' @param params <[`tidy-select`][tidyselect::language]> parameters for the
#'   trace and density plots (default: the first `max_params`).
#' @param max_params number of parameters shown when `params` is `NULL`.
#' @param ... unused.
#' @return a `ggplot` object.
#' @family Bayesian VAR
#' @importFrom rlang .data
#' @examples
#' set.seed(1)
#' d <- mcmc_diagnostics(cbind(a = rnorm(500), b = rnorm(500)))
#' autoplot(d)
#' autoplot(d, "trace")
#' @exportS3Method ggplot2::autoplot
autoplot.mcmc_diagnostics <- function(
  object,
  type = c("summary", "trace", "density"),
  params = NULL,
  max_params = 12L,
  ...
) {
  type <- match.arg(type)
  if (type == "summary") {
    return(diagnostics_summary_plot(tidy.mcmc_diagnostics(object)))
  }
  a <- attr(object, "draws")
  nm <- dimnames(a)[[3]]
  keep <- select_or(
    rlang::enquo(params),
    nm,
    "params",
    seq_len(min(max_params, length(nm)))
  )
  d <- chain_long(a[,, keep, drop = FALSE])
  if (type == "trace") {
    return(trace_plot(d))
  }
  ggplot2::ggplot(d, ggplot2::aes(.data$value, colour = .data$chain)) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap("parameter", scales = "free") +
    ggplot2::scale_colour_manual(
      values = varm_palette(dim(a)[2]),
      name = "Chain"
    ) +
    ggplot2::labs(x = NULL, y = "Posterior density") +
    theme_varm()
}

#' @rdname autoplot.mcmc_diagnostics
#' @param x an `mcmc_diagnostics` object.
#' @export
plot.mcmc_diagnostics <- function(x, ...) {
  p <- ggplot2::autoplot(x, ...)
  print(p)
  invisible(p)
}

# long tibble (iteration, chain, parameter, value) of an iter x chain x param array
chain_long <- function(a) {
  d <- dim(a)
  tibble::tibble(
    iteration = rep(seq_len(d[1]), d[2] * d[3]),
    chain = factor(rep(rep(seq_len(d[2]), each = d[1]), d[3])),
    parameter = factor(
      rep(dimnames(a)[[3]], each = d[1] * d[2]),
      levels = dimnames(a)[[3]]
    ),
    value = as.vector(a)
  )
}

trace_plot <- function(d) {
  ggplot2::ggplot(
    d,
    ggplot2::aes(.data$iteration, .data$value, colour = .data$chain)
  ) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.8) +
    ggplot2::facet_wrap("parameter", scales = "free_y") +
    ggplot2::scale_colour_manual(
      values = varm_palette(nlevels(d$chain)),
      name = "Chain"
    ) +
    ggplot2::labs(x = "Iteration", y = NULL) +
    theme_varm()
}

diagnostics_summary_plot <- function(d) {
  long <- tibble::tibble(
    term = factor(rep(d$term, 2L), levels = rev(d$term)),
    statistic = rep(c("Effective sample size", "R-hat"), each = nrow(d)),
    value = c(d$ess, d$rhat)
  )
  ggplot2::ggplot(long, ggplot2::aes(.data$value, .data$term)) +
    ggplot2::geom_point(colour = varm_palette(1), size = 2) +
    ggplot2::facet_wrap("statistic", scales = "free_x") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_varm()
}

as_chain_array <- function(d) {
  if (inherits(d, "bvarm")) {
    d <- if (!is.null(d$hyper)) d$hyper else coef_draw_matrix(d$draws$B)
  } else if (inherits(d, "varm_draws")) {
    d <- coef_draw_matrix(d$B)
  }
  if (is.list(d) && !is.data.frame(d)) {
    return(chains_from_list(d))
  }
  if (is.data.frame(d)) {
    d <- as.matrix(d)
  }
  if (is.null(dim(d))) {
    d <- matrix(d, ncol = 1)
  }
  if (length(dim(d)) == 2L) {
    return(array(
      d,
      c(nrow(d), 1L, ncol(d)),
      dimnames = list(NULL, NULL, colnames(d))
    ))
  }
  if (length(dim(d)) == 3L) {
    return(d)
  }
  rlang::abort(
    "`draws` must be a matrix, a 3-d array (iter x chain x param), a list of matrices, or a bvarm/varm_draws"
  )
}

# iterations x parameters matrix from a coefficient array (m x K x N)
coef_draw_matrix <- function(B) {
  nm <- NULL
  if (!is.null(rownames(B)) && !is.null(colnames(B))) {
    nm <- as.vector(outer(rownames(B), colnames(B), paste, sep = ":"))
  }
  out <- t(matrix(B, dim(B)[1] * dim(B)[2]))
  colnames(out) <- nm
  out
}

chains_from_list <- function(d) {
  d <- lapply(d, as.matrix)
  n <- min(vapply(d, nrow, 0L))
  a <- array(
    NA_real_,
    c(n, length(d), ncol(d[[1]])),
    dimnames = list(NULL, NULL, colnames(d[[1]]))
  )
  for (j in seq_along(d)) {
    a[, j, ] <- d[[j]][seq_len(n), , drop = FALSE]
  }
  a
}

# autocovariance (biased, divisor n) of one chain via FFT
acov_fft <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  nf <- 2L^ceiling(log2(2 * n))
  f <- stats::fft(c(x, rep(0, nf - n)))
  Re(stats::fft(f * Conj(f), inverse = TRUE))[seq_len(n)] / nf / n
}

# Vehtari et al. (2021) multi-chain ESS with Geyer's initial monotone positive
# sequence
ess_chains <- function(x) {
  n <- nrow(x)
  m <- ncol(x)
  if (stats::var(as.vector(x)) == 0) {
    return(NA_real_)
  }
  g <- matrix(
    vapply(seq_len(m), function(j) acov_fft(x[, j]), numeric(n)),
    n,
    m
  )
  W <- mean(apply(x, 2, stats::var))
  Bn <- if (m > 1) stats::var(colMeans(x)) else 0 # B/n
  vplus <- W * (n - 1) / n + Bn
  if (vplus <= 0) {
    return(NA_real_)
  }
  gm <- rowMeans(g) * n / (n - 1) # unbiased-ish scaling consistent with W
  rho <- 1 - (W - gm) / vplus
  rho[1] <- 1
  npair <- (n - 1) %/% 2
  # Geyer pairs (lags 2t, 2t+1)
  pr <- rho[2 * seq_len(npair) - 1] + rho[2 * seq_len(npair)]
  neg <- which(pr <= 0)
  if (length(neg)) {
    pr <- pr[seq_len(neg[1] - 1)]
  }
  if (length(pr) > 1) {
    pr <- cummin(pr) # monotone
  }
  tau <- -1 + 2 * sum(pr)
  max(1, m * n / max(tau, 1 / log10(m * n)))
}

rhat_split <- function(x) {
  n <- nrow(x)
  h <- n %/% 2
  s <- cbind(
    x[seq_len(h), , drop = FALSE],
    x[n - h + seq_len(h), , drop = FALSE]
  )
  W <- mean(apply(s, 2, stats::var))
  Bn <- stats::var(colMeans(s)) # B/h
  if (W == 0) {
    return(NA_real_)
  }
  sqrt(((h - 1) / h * W + Bn) / W)
}

geweke_z <- function(x, f1 = 0.1, f2 = 0.5) {
  n <- length(x)
  a <- x[seq_len(floor(f1 * n))]
  b <- x[(n - floor(f2 * n) + 1):n]
  va <- stats::var(a) / ess_chains(matrix(a))
  vb <- stats::var(b) / ess_chains(matrix(b))
  (mean(a) - mean(b)) / sqrt(va + vb)
}
