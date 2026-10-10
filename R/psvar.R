#' Proxy-SVAR (external instrument) identification
#'
#' Mertens & Ravn (2013) / Stock & Watson (2012) proxy identification of one structural
#' shock in a reduced-form [varm()] fit. The impact vector is \eqn{\gamma = T^{-1}\sum_t u_t m_t}
#' (covariance of the reduced-form residuals with the proxy), scaled so that the impact
#' response of `target` equals `size`. Responses are \eqn{\Theta_h \gamma}, reusing the
#' MA machinery of [irf()]. Bands come from the proxy wild bootstrap (Rademacher weights
#' applied jointly to the residuals and the proxy, Jentsch & Lunsford 2019).
#'
#' This does not fit the `id_*` interface because [irf()] maps `(B, Sigma)` to impact
#' matrices, whereas a proxy also needs the residuals and the instrument.
#'
#' @param fit a [varm()] object.
#' @param m proxy: numeric vector (or one-column matrix) of length `nrow(fit$Y)`, aligned with
#'   the data (the first `p` values are dropped to match the residuals).
#' @param target variable (name or position) whose impact response is normalised to `size`.
#' @param size size of the shock (impact response of `target`).
#' @param horizon maximum horizon (impact is 0).
#' @param nboot number of wild-bootstrap draws (`0` for point estimates only).
#' @param level confidence levels; bands are equal-tailed percentiles.
#' @return a `varm_proxy` object: `estimate` (`(H+1) x K`), `draws` (`(H+1) x K x nboot`),
#'   `bands` (`(H+1) x K x 2*length(level)`), `gamma`, `fstat` (see [proxy_fstat()]).
#' @examples
#' set.seed(1)
#' e <- matrix(rnorm(400), 200)
#' y <- matrix(0, 200, 2)
#' for (t in 2:200) y[t, ] <- 0.5 * y[t - 1, ] + c(e[t, 1] + 0.3 * e[t, 2], e[t, 2] + 0.5 * e[t, 1])
#' fit <- varm(y, p = 1)
#' varm_proxy(fit, m = e[, 1] + rnorm(200, sd = 0.5), horizon = 8, nboot = 100)
#' @export
varm_proxy <- function(
  fit,
  m,
  target = 1L,
  size = 1,
  horizon = 20L,
  nboot = 500L,
  level = c(0.68, 0.90)
) {
  if (!inherits(fit, "varm")) {
    stop("`fit` must be a varm object", call. = FALSE)
  }
  p <- fit$p
  K <- fit$K
  Tn <- nrow(fit$Y)
  m <- as.numeric(m)
  if (length(m) != Tn || anyNA(m)) {
    stop("`m` must have no NAs and length nrow(fit$Y)", call. = FALSE)
  }
  j <- resolve_var(target, fit$varnames, "target")
  me <- m[-seq_len(p)]
  proxy_irf <- function(B, U, me) {
    g <- drop(crossprod(U, me)) / nrow(U)
    g <- size * g / g[j]
    Th <- cpp_ma(B, K, p, horizon)
    t(vapply(
      seq_len(horizon + 1L),
      function(h) drop(Th[,, h] %*% g),
      numeric(K)
    ))
  }
  est <- proxy_irf(fit$B, fit$U, me)
  g <- drop(crossprod(fit$U, me)) / nrow(fit$U)
  dimnames(est) <- list(0:horizon, fit$varnames)
  out <- list(
    estimate = est,
    gamma = setNames(g, fit$varnames),
    horizon = horizon,
    p = p,
    K = K,
    varnames = fit$varnames,
    target = j,
    size = size,
    level = level,
    fstat = proxy_fstat(fit, m, j)
  )
  if (nboot > 0L) {
    Tp <- nrow(fit$U)
    detp <- fit$det[-seq_len(p), , drop = FALSE]
    Y0 <- fit$Y[seq_len(p), , drop = FALSE]
    bd <- array(0, c(horizon + 1L, K, nboot))
    for (b in seq_len(nboot)) {
      wb <- sample(c(-1, 1), Tp, replace = TRUE)
      Yb <- rbind(Y0, cpp_path(fit$B, K, p, Y0, detp, fit$U * wb))
      e <- cpp_var_ols(Yb, p, fit$det)
      bd[,, b] <- proxy_irf(e$B, e$U, me * wb)
    }
    probs <- level_probs(level)
    out$draws <- bd
    out$bands <- aperm(
      apply(bd, c(1, 2), stats::quantile, probs = probs),
      c(2, 3, 1)
    )
    dimnames(out$bands) <- list(0:horizon, fit$varnames, format(probs))
  }
  structure(out, class = "varm_proxy")
}

#' First-stage strength of a proxy
#'
#' Three F statistics for the relevance of the proxy for the `target` shock:
#' `F` from regressing the reduced-form residual of `target` on the proxy,
#' `F_robust` its heteroskedasticity-robust (HC1) version (needs \pkg{sandwich},
#' else `NA`), and `F_hac`, the Montiel Olea-Stock-Watson (2021) Wald statistic
#' \eqn{T\hat\Gamma_j^2/\hat W_{\Gamma_{jj}}} using the delta-method HAC covariance
#' that also accounts for the estimation of the VAR coefficients. `F_hac` is the
#' relevant statistic for [varm_proxy_robust()].
#' @inheritParams varm_proxy
#' @param nwlags Newey-West truncation lag for `F_hac`; default `floor(4 * (T/100)^(2/9))`.
#' @return named numeric `c(F = , F_robust = , F_hac = )`.
#' @seealso [varm_proxy_robust()]
#' @export
proxy_fstat <- function(fit, m, target = 1L, nwlags = NULL) {
  j <- resolve_var(target, fit$varnames, "target")
  me <- as.numeric(m)[-seq_len(fit$p)]
  mod <- stats::lm(fit$U[, j] ~ me)
  f <- unname(summary(mod)$fstatistic[1])
  fr <- if (requireNamespace("sandwich", quietly = TRUE)) {
    unname((stats::coef(mod)[2]^2) / sandwich::vcovHC(mod, type = "HC1")[2, 2])
  } else {
    NA_real_
  }
  cv <- proxy_cov(fit, me, nwlags)
  fh <- cv$Tp * cv$Gamma[j]^2 / cv$W[cv$iG[j], cv$iG[j]]
  c(F = f, F_robust = fr, F_hac = unname(fh))
}

#' @export
print.varm_proxy <- function(x, ...) {
  cat(
    "<varm_proxy> shock to ",
    x$varnames[x$target],
    " (size ",
    x$size,
    "), horizon 0..",
    x$horizon,
    if (!is.null(x$draws)) {
      paste0(", ", dim(x$draws)[3], " wild-bootstrap draws")
    },
    "\n",
    sep = ""
  )
  cat(
    "first-stage F = ",
    format(x$fstat[["F"]], digits = 4),
    "; robust F = ",
    format(x$fstat[["F_robust"]], digits = 4),
    "\n",
    sep = ""
  )
  invisible(x)
}
