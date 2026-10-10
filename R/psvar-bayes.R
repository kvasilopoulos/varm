#' Bayesian proxy-SVAR on top of [bvarm()] draws
#'
#' A simple, limited-information Bayesian proxy-SVAR. For each posterior draw
#' \eqn{(B_s, \Sigma_s)} of a [bvarm()] fit, with residuals \eqn{u_t^{(s)}}, the proxy is
#' modelled by the measurement regression \eqn{m_t = \Lambda' u_t + e_t},
#' \eqn{e_t \sim N(0, \Sigma_e)}, with a flat (Jeffreys) prior on \eqn{(\Lambda, \Sigma_e)}.
#' A draw of \eqn{(\Lambda_s, \Sigma_{e,s})} from its (matrix-normal / inverse-Wishart)
#' conditional posterior gives \eqn{\mathrm{Cov}(u, m) = \Sigma_s\Lambda_s} and
#' \eqn{\mathrm{Var}(m) = \Lambda_s'\Sigma_s\Lambda_s + \Sigma_{e,s}}, from which the
#' impact vector is obtained exactly as in [varm_proxy()] (also for `r > 1`) and the
#' responses are \eqn{\Theta_h(B_s)\,\gamma_s}.
#'
#' **Limitations (read before using).** This is *not* the joint-likelihood Bayesian
#' proxy-SVAR of Caldara & Herbst (2019) or Arias, Rubio-Ramirez & Waggoner (2021): the
#' proxy does not enter the likelihood of the VAR, so it does not update \eqn{(B,\Sigma)}
#' and the reduced-form posterior is whatever the [bvarm()] prior implies; there is no prior
#' on the structural parameters or on the instrument-shock relationship other than the flat
#' prior on \eqn{\Lambda}; stochastic-volatility fits (time-varying covariance) are rejected. With a weak instrument the posterior
#' of the normalised responses is correctly wide (heavy-tailed) but, being Bayesian with a
#' flat prior, its credible sets do not have frequentist coverage guarantees
#' (use [varm_proxy_robust()] for those). Credible bands are equal-tailed posterior
#' percentiles; `estimate` is the posterior median.
#' @inheritParams varm_proxy
#' @param fit a [bvarm()] object (class `bvarm`).
#' @param ndraw number of posterior draws to use (default: all stored draws).
#' @return a `varm_proxy` object with `method = "bayes"` (same layout as [varm_proxy()],
#'   `draws` are posterior draws).
#' @references Caldara, D. and Herbst, E. (2019). Monetary policy, real activity, and credit
#'   spreads: evidence from Bayesian proxy SVARs. *AEJ: Macroeconomics* 11(1), 157-192.
#'
#'   Arias, J. E., Rubio-Ramirez, J. F. and Waggoner, D. F. (2021). Inference in Bayesian
#'   proxy-SVARs. *Journal of Econometrics* 225(1), 88-106.
#'
#'   Giacomini, R., Kitagawa, T. and Read, M. (2022). Robust Bayesian inference in proxy
#'   SVARs. *Journal of Econometrics* (doi:10.1016/j.jeconom.2021.02.003); not implemented.
#' @seealso [varm_proxy()], [varm_proxy_robust()]
#' @examples
#' \donttest{
#' set.seed(1)
#' e <- matrix(rnorm(400), 200)
#' y <- matrix(0, 200, 2)
#' for (t in 2:200) y[t, ] <- 0.5 * y[t - 1, ] + c(e[t, 1] + 0.3 * e[t, 2], e[t, 2] + 0.5 * e[t, 1])
#' bf <- bvarm(y, p = 1)
#' varm_proxy_bayes(bf, m = e[, 1] + rnorm(200, sd = 0.5), horizon = 8, ndraw = 200)
#' }
#' @export
varm_proxy_bayes <- function(
  fit,
  m,
  target = 1L,
  size = 1,
  horizon = 20L,
  level = c(0.68, 0.90),
  ndraw = NULL
) {
  if (!inherits(fit, "bvarm")) {
    stop("`fit` must be a bvarm object", call. = FALSE)
  }
  d <- fit$draws
  Bs <- d$B
  Ss <- d$Sigma
  if (length(dim(Ss)) != 3L || !is.null(fit$Sigma_mean)) {
    stop(
      "`fit` has a constant innovation covariance (stochastic-volatility fits are not supported)",
      call. = FALSE
    )
  }
  p <- fit$p
  K <- fit$K
  H <- as.integer(horizon)
  M <- proxy_check_m(m, nrow(fit$Y))[-seq_len(p), , drop = FALSE]
  r <- ncol(M)
  iP <- resolve_var(target, fit$varnames, "target")
  if (length(iP) != r || anyDuplicated(iP)) {
    stop(
      "`target` must have one distinct variable per instrument",
      call. = FALSE
    )
  }
  size <- rep_len(size, r)
  Tp <- nrow(M)
  if (Tp - K <= 0L) {
    stop("too few observations for the measurement regression", call. = FALSE)
  }
  N <- dim(Bs)[3]
  idx <- if (is.null(ndraw) || ndraw >= N) {
    seq_len(N)
  } else {
    sort(sample.int(N, ndraw))
  }
  X <- proxy_X(fit)
  Yt <- fit$Y[-seq_len(p), , drop = FALSE]
  bd <- array(0, c(H + 1L, K, r, length(idx)))
  for (k in seq_along(idx)) {
    s <- idx[k]
    U <- Yt - X %*% Bs[,, s]
    UUi <- solve(crossprod(U))
    Lh <- UUi %*% crossprod(U, M)
    E <- M - U %*% Lh
    # Sigma_e ~ IW(Tp - K, E'E), Lambda | Sigma_e ~ MN(Lh, (U'U)^-1, Sigma_e)
    Se <- solve(stats::rWishart(1L, Tp - K, solve(crossprod(E)))[,, 1])
    Lam <- Lh + t(chol(UUi)) %*% matrix(stats::rnorm(K * r), K, r) %*% chol(Se)
    Sg <- Ss[,, s]
    mom <- list(
      Su = Sg,
      Sum = Sg %*% Lam,
      Smm = crossprod(Lam, Sg %*% Lam) + Se
    )
    bd[,,, k] <- proxy_irf_mom(Bs[,, s], mom, iP, size, K, p, H)
  }
  sn <- if (r == 1L) NULL else paste0("shock", seq_len(r))
  out <- list(
    horizon = H,
    p = p,
    K = K,
    r = r,
    varnames = fit$varnames,
    target = iP,
    size = size,
    level = level,
    method = "bayes"
  )
  out <- proxy_bands(out, bd, level, fit$varnames, sn)
  out$estimate <- proxy_drop(
    apply(bd, 1:3, stats::median),
    r,
    list(0:H, fit$varnames, sn)
  )
  structure(out, class = "varm_proxy")
}
