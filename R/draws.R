#' Reduced-form parameter draws
#'
#' Common currency for uncertainty: a `varm_draws` object holds `N` draws of the
#' coefficient matrix `B` (m x K x N) and innovation covariance `Sigma` (K x K x N).
#' Bootstrap, asymptotic Monte Carlo and every Bayesian sampler in \pkg{varm} produce one,
#' so identification, IRFs, FEVDs and forecasts share one code path.
#'
#' @param object a `varm` fit.
#' @param method `"bootstrap"` (recursive-design residual bootstrap), `"wild"`
#'   (Rademacher wild bootstrap, robust to conditional heteroskedasticity),
#'   `"wild_mammen"` or `"mc"` (asymptotic normal / inverse-Wishart Monte Carlo).
#' @param n number of draws.
#' @param bias_correct apply Kilian (1998) bootstrap-after-bootstrap bias correction.
#' @export
draws <- function(object, ...) UseMethod("draws")

#' @rdname draws
#' @export
draws.varm <- function(
  object,
  method = c("bootstrap", "wild", "wild_mammen", "mc"),
  n = 1000L,
  bias_correct = FALSE,
  ...
) {
  method <- match.arg(method)
  n <- as.integer(n)
  K <- object$K
  p <- object$p
  if (method == "mc") {
    d <- cpp_var_mc(object$B, object$U, object$XtXinv, n)
  } else {
    type <- c(bootstrap = 0L, wild = 1L, wild_mammen = 2L)[[method]]
    B0 <- object$B
    U0 <- object$U
    bias <- NULL
    if (bias_correct) {
      d1 <- cpp_var_boot(object$Y, p, object$det, B0, U0, n, type)
      bias <- apply(d1$B, c(1, 2), mean) - B0
      B0 <- stationary_adjust(B0 - bias, B0, K, p, bias)
      U0 <- object$Y[-seq_len(p), , drop = FALSE] - cpp_var_ols_X(object, B0)
    }
    d <- cpp_var_boot(object$Y, p, object$det, B0, U0, n, type)
    if (bias_correct) {
      for (i in seq_len(n)) {
        d$B[,, i] <- stationary_adjust(d$B[,, i] - bias, d$B[,, i], K, p, bias)
      }
    }
  }
  structure(
    list(
      B = d$B,
      Sigma = d$Sigma,
      p = p,
      K = K,
      varnames = object$varnames,
      method = method,
      fit = object
    ),
    class = "varm_draws"
  )
}

cpp_var_ols_X <- function(object, B) {
  X <- cpp_var_ols(object$Y, object$p, object$det)$X
  X %*% B
}

# Kilian (1998) stationarity correction: shrink the bias until the VAR is stable
stationary_adjust <- function(Bc, B, K, p, bias) {
  if (max(cpp_stability(Bc, K, p)) < 1) {
    return(Bc)
  }
  delta <- 1
  while (delta > 0 && max(cpp_stability(Bc, K, p)) >= 1) {
    delta <- delta - 0.01
    Bc <- B - delta * bias
  }
  Bc
}

#' @export
print.varm_draws <- function(x, ...) {
  cat(
    "<varm_draws> ",
    dim(x$B)[3],
    " draws (",
    x$method,
    "); K = ",
    x$K,
    ", p = ",
    x$p,
    "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
stability.varm_draws <- function(object) {
  t(apply(object$B, 3, function(b) cpp_stability(b, object$K, object$p)))
}
