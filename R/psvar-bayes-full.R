# Full-posterior Bayesian proxy-SVARs: ARW-type exact draws and GKR robust Bayes.
# Shared machinery: first-stage regression posterior + rotation sampler cpp_proxy_rot().

# ---- shared internals --------------------------------------------------------

# Posterior ingredients of the first-stage regression (GKR 2022, eq. 10; ARW 2021 lags-of-y case)
#   m_t = D y_t + G x_t + sum_l H_l m_{t-l} + v_t,   Var(v_t) = Upsilon
# with Jeffreys prior |Upsilon|^{-(r+1)/2}, flat on the coefficients.
proxy_first_stage <- function(fit, m, pm = 0L) {
  p <- fit$p
  K <- fit$K
  M <- proxy_check_m(m, nrow(fit$Y))[-seq_len(p), , drop = FALSE]
  Te <- nrow(M)
  r <- ncol(M)
  pm <- as.integer(pm)
  Z <- cbind(fit$Y[-seq_len(p), , drop = FALSE], proxy_X(fit))
  if (pm > 0L) {
    Z <- cbind(
      Z,
      do.call(
        cbind,
        lapply(seq_len(pm), function(l) {
          M[(pm + 1L - l):(Te - l), , drop = FALSE]
        })
      )
    )
  }
  rows <- (pm + 1L):Te
  Z <- Z[rows, , drop = FALSE]
  Mr <- M[rows, , drop = FALSE]
  df <- nrow(Z) - ncol(Z)
  if (df <= r) {
    stop("too few observations for the first-stage regression", call. = FALSE)
  }
  ZZi <- solve(crossprod(Z))
  Ch <- ZZi %*% crossprod(Z, Mr)
  E <- Mr - Z %*% Ch
  list(K = K, r = r, Ch = Ch, L = t(chol(ZZi)), S = crossprod(E), df = df)
}

# N draws of D (r x K x N) and Upsilon (r x r x N)
proxy_first_stage_draw <- function(fs, N) {
  r <- fs$r
  q <- nrow(fs$Ch)
  W <- stats::rWishart(N, fs$df, solve(fs$S))
  D <- array(0, c(r, fs$K, N))
  Ups <- array(0, c(r, r, N))
  for (n in seq_len(N)) {
    U <- solve(matrix(W[,, n], r, r))
    Cf <- fs$Ch + fs$L %*% matrix(stats::rnorm(q * r), q, r) %*% chol(U)
    D[,, n] <- t(Cf[seq_len(fs$K), , drop = FALSE])
    Ups[,, n] <- U
  }
  list(D = D, Ups = Ups)
}

# K x K x Hs sign cube; `S` is K x ncol(cols)(xH) for the proxy shocks (cols = last r columns) or K x K(xH)
proxy_sign_cube <- function(S, K, r, horizon, proxy_only) {
  if (is.null(S)) {
    return(array(0, c(K, K, 1L)))
  }
  nc <- if (proxy_only) r else K
  if (length(dim(S)) == 2L) {
    S <- array(rep(S, horizon), c(dim(S), horizon))
  }
  if (length(dim(S)) != 3L || dim(S)[1] != K || dim(S)[2] != nc) {
    stop(
      "`signs` must be K x ",
      nc,
      " or K x ",
      nc,
      " x H",
      if (proxy_only) " (columns = the r proxy shocks)",
      call. = FALSE
    )
  }
  S[is.na(S)] <- 0
  if (!all(S %in% c(-1, 0, 1))) {
    stop("`signs` entries must be -1, 0 or +1", call. = FALSE)
  }
  if (!proxy_only) {
    return(S)
  }
  out <- array(0, c(K, K, dim(S)[3]))
  out[, K - r + seq_len(r), ] <- S
  out
}

# minimum eigenvalue of the reliability matrix (Upsilon + D Sigma D')^{-1} D Sigma D'  (ARW 2021, sec. 2.3)
proxy_reliability <- function(D, Ups, Sigma) {
  V <- D %*% Sigma %*% t(D)
  min(Re(eigen(solve(Ups + V, V), only.values = TRUE)$values))
}

proxy_draws_check <- function(fit) {
  if (!inherits(fit, "bvarm")) {
    stop("`fit` must be a bvarm object", call. = FALSE)
  }
  if (length(dim(fit$draws$Sigma)) != 3L || !is.null(fit$Sigma_mean)) {
    stop(
      "`fit` has a constant innovation covariance (stochastic-volatility fits are not supported)",
      call. = FALSE
    )
  }
}

# ---- (A) ARW-type posterior --------------------------------------------------

#' Bayesian proxy-SVAR with exogeneity, relevance and sign restrictions (Arias-Rubio-Ramirez-Waggoner type)
#'
#' Full-likelihood Bayesian proxy-SVAR in the spirit of Arias, Rubio-Ramirez & Waggoner (2021).
#' The `r` proxies obey the exogeneity restrictions \eqn{E[m_t\varepsilon_{jt}'] = 0} for the first
#' \eqn{K-r} shocks and a relevance condition for the **last `r`** shocks. As in their
#' Proxy-SVAR (their eq. 3 with the block restrictions), the data have the joint likelihood
#' \eqn{p(Y)\,p(M\mid Y)} with the proxy equation
#' \deqn{m_t = D y_t + G x_t + \sum_{l=1}^{p_m} H_l m_{t-l} + v_t,\quad v_t\sim N(0,\Upsilon),}
#' which the lags of \eqn{y_t} (and optionally of \eqn{m_t}) are allowed to enter. The likelihood
#' factorises, so a posterior draw is: \eqn{(B,\Sigma)} from the supplied [bvarm()] fit (any
#' conjugate prior, e.g. Minnesota) and \eqn{(D,G,H,\Upsilon)} from the first-stage regression
#' (Jeffreys prior), independently. The exogeneity restrictions
#' \eqn{D\Sigma_{tr}q_j=0,\ j\le K-r}, with \eqn{\Sigma=\Sigma_{tr}\Sigma_{tr}'} and \eqn{Q}
#' the rotation matrix, are linear in the columns of \eqn{Q}; \eqn{Q} is drawn column by column
#' from the unit sphere of the admissible subspace (ARW Algorithm 1 / Arias et al. 2018, using
#' the orthogonal-reduced-form parameterisation), and sign restrictions on the impulse responses of
#' the proxy shocks are imposed by accept/reject. With `r = 1` the proxy shock is
#' **point-identified** (\eqn{q_K \propto \Sigma_{tr}'D'}) and signs only truncate the posterior; with
#' `r > 1` exogeneity alone leaves the proxy shocks identified only up to an orthogonal
#' \eqn{r\times r} rotation, which the sign restrictions (needed for identification) pin down.
#' A reliability truncation \eqn{\lambda_{\min}\{(\Upsilon+D\Sigma D')^{-1}D\Sigma D'\}\ge\gamma}
#' (ARW 2021 section 2.3, after Caldara & Herbst 2019) is available through `reliability`.
#'
#' **What is and is not the ARW posterior.** The prior implied is the conjugate prior on the
#' reduced form times the (conditionally) uniform prior on the admissible rotations, with
#' importance weights set to one. ARW's Algorithm 2 instead targets a normal-generalised-normal
#' prior on the *structural* parameters and re-weights draws by a volume element
#' (their Theorem 3, taken from Arias et al. 2018); that element is not given in the paper and is **not
#' implemented here**. The two posteriors coincide when the proxy shocks are point-identified
#' (`r = 1`, no sign restrictions on non-proxy shocks) up to the reduced-form prior; with
#' `r > 1` they differ through the prior on the rotation. Non-proxy shocks are not identified
#' and their IRF columns are returned as zero. Compare [varm_proxy_bayes()], which
#' does not allow lags of \eqn{y_t} in the proxy equation and has no sign or reliability
#' restrictions, and [varm_proxy_gkr()] for prior-robust inference when `r > 1` or shocks other than
#' the proxy shocks are of interest.
#'
#' Not implemented: Caldara & Herbst (2019)'s Metropolis-within-Gibbs sampler with a
#' structural measurement equation \eqn{m_t=\beta\varepsilon_{1t}+\sigma_\nu\nu_t} (where the
#' proxy likelihood depends on \eqn{(B,\Sigma)} through \eqn{u_t}, so the proxy updates the VAR
#' coefficients) and priors on \eqn{\beta,\sigma_\nu}.
#'
#' @inheritParams varm_proxy_bayes
#' @param signs sign restrictions on the responses of the proxy shocks: `K x r` matrix (or
#'   `K x r x H` array) of `+1`, `-1`, `0`; column `i` is the `i`-th proxy shock (shock `K - r + i`).
#'   `NULL` for none. The sign of a shock is otherwise normalised by `size` (impact on `target`).
#' @param size impact response of `target` to each proxy shock (recycled), or `NULL` for
#'   one-standard-deviation shocks (sign normalised so that \eqn{A_0} has non-negative diagonal;
#'   state the sign through `signs` in that case).
#' @param reliability minimum reliability \eqn{\gamma\in[0,1)}; draws below are discarded (truncated prior).
#' @param sign_horizon number of horizons (from impact) over which a *matrix* `signs` is imposed.
#' @param pm number of lags of the proxy in the first-stage regression.
#' @param nrot accepted rotations kept per reduced-form draw.
#' @param maxtry maximum rotation draws per reduced-form draw.
#' @return `varm_proxy_arw()`: a `varm_proxy` object (`method = "arw"`) with `estimate` (posterior median),
#'   `draws`, `bands`, plus `accept_rate`, `reliability` (posterior median of the reliability
#'   index, per proxy). `id_proxy_bayes()`: a `varm_id` for [irf()] on a [bvarm()] fit (columns of
#'   non-proxy shocks are zero).
#' @references Arias, J. E., Rubio-Ramirez, J. F. and Waggoner, D. F. (2021). Inference in Bayesian
#'   proxy-SVARs. *Journal of Econometrics* 225(1), 88-106.
#'
#'   Arias, J. E., Rubio-Ramirez, J. F. and Waggoner, D. F. (2018). Inference based on SVARs
#'   identified with sign and zero restrictions. *Econometrica* 86(2), 685-720.
#'
#'   Caldara, D. and Herbst, E. (2019). Monetary policy, real activity, and credit spreads:
#'   evidence from Bayesian proxy SVARs. *AEJ: Macroeconomics* 11(1), 157-192.
#' @seealso [varm_proxy_gkr()], [varm_proxy_bayes()], [varm_proxy()]
#' @examples
#' \donttest{
#' set.seed(1)
#' e <- matrix(rnorm(900), 300)
#' y <- matrix(0, 300, 3)
#' for (t in 2:300) y[t, ] <- 0.4 * y[t - 1, ] + e[t, ] %*% t(matrix(c(1, .3, .2, 0, 1, .4, .5, .2, 1), 3))
#' bf <- bvarm(y, p = 1, ndraw = 300, nburn = 100)
#' varm_proxy_arw(bf, m = e[, 3] + rnorm(300, sd = 0.5), target = 3, horizon = 8)
#' }
#' @export
varm_proxy_arw <- function(
  fit,
  m,
  target = 1L,
  size = 1,
  signs = NULL,
  reliability = NULL,
  pm = 0L,
  horizon = 20L,
  level = c(0.68, 0.90),
  ndraw = NULL,
  nrot = 1L,
  maxtry = 1000L,
  sign_horizon = 1L
) {
  proxy_draws_check(fit)
  K <- fit$K
  p <- fit$p
  H <- as.integer(horizon)
  id <- id_proxy_bayes(
    fit,
    m,
    target,
    size,
    signs,
    reliability,
    pm,
    nrot,
    maxtry,
    sign_horizon
  )
  r <- id$r
  Bs <- fit$draws$B
  Ss <- fit$draws$Sigma
  N <- dim(Bs)[3]
  idx <- if (is.null(ndraw) || ndraw >= N) {
    seq_len(N)
  } else {
    sort(sample.int(N, ndraw))
  }
  si <- id$structural(
    Bs[,, idx, drop = FALSE],
    Ss[,, idx, drop = FALSE],
    K,
    p,
    H,
    id,
    fit$varnames
  )
  irf4 <- si$irf[,,,, drop = FALSE]
  Na <- dim(irf4)[4]
  bd <- array(0, c(H + 1L, K, r, Na))
  for (j in seq_len(r)) {
    bd[,, j, ] <- aperm(irf4[, K - r + j, , ], c(2, 1, 3))
  }
  iP <- resolve_var(target, fit$varnames, "target")
  sn <- if (r == 1L) NULL else paste0("shock", seq_len(r))
  out <- list(
    horizon = H,
    p = p,
    K = K,
    r = r,
    varnames = fit$varnames,
    target = iP,
    size = if (is.null(size)) NA_real_ else rep_len(size, r),
    level = level,
    method = "arw",
    accept_rate = si$accept_rate,
    nacc = Na
  )
  out <- proxy_bands(out, bd, level, fit$varnames, sn)
  out$estimate <- proxy_drop(
    apply(bd, 1:3, stats::median),
    r,
    list(0:H, fit$varnames, sn)
  )
  structure(out, class = "varm_proxy")
}

#' @rdname varm_proxy_arw
#' @export
id_proxy_bayes <- function(
  fit,
  m,
  target = 1L,
  size = 1,
  signs = NULL,
  reliability = NULL,
  pm = 0L,
  nrot = 1L,
  maxtry = 1000L,
  sign_horizon = 1L
) {
  fs <- proxy_first_stage(fit, m, pm)
  r <- fs$r
  K <- fit$K
  iP <- resolve_var(target, fit$varnames, "target")
  if (length(iP) != r || anyDuplicated(iP)) {
    stop(
      "`target` must have one distinct variable per instrument",
      call. = FALSE
    )
  }
  if (r >= K) {
    stop("need fewer instruments than variables", call. = FALSE)
  }
  if (
    !is.null(reliability) &&
      (length(reliability) != 1L || reliability < 0 || reliability >= 1)
  ) {
    stop("`reliability` must be a number in [0, 1)", call. = FALSE)
  }
  structure(
    list(
      type = "proxy_bayes",
      fs = fs,
      r = r,
      target = iP,
      size = if (is.null(size)) NULL else rep_len(size, r),
      S = proxy_sign_cube(signs, K, r, sign_horizon, TRUE),
      reliability = reliability,
      nrot = as.integer(nrot),
      maxtry = as.integer(maxtry),
      structural = proxy_bayes_structural
    ),
    class = c("id_proxy_bayes", "varm_id")
  )
}

proxy_bayes_structural <- function(Bs, Ss, K, p, H, id, varnames) {
  N <- dim(Bs)[3]
  r <- id$r
  fd <- proxy_first_stage_draw(id$fs, N)
  keep <- seq_len(N)
  rel <- NULL
  if (!is.null(id$reliability)) {
    rel <- vapply(
      keep,
      function(n) {
        proxy_reliability(
          matrix(fd$D[,, n], r, K),
          matrix(fd$Ups[,, n], r, r),
          matrix(Ss[,, n], K, K)
        )
      },
      0
    )
    keep <- which(rel >= id$reliability)
    if (!length(keep)) {
      stop("no draw satisfied the reliability threshold", call. = FALSE)
    }
  }
  o <- cpp_proxy_rot(
    Bs[,, keep, drop = FALSE],
    Ss[,, keep, drop = FALSE],
    fd$D[,, keep, drop = FALSE],
    p,
    H,
    id$S,
    id$nrot,
    id$nrot,
    id$maxtry,
    -1L,
    -1L,
    0
  )
  if (!length(o$draw)) {
    stop(
      "no rotation satisfied the sign restrictions (tries = ",
      o$tries,
      ")",
      call. = FALSE
    )
  }
  H1 <- H + 1L
  cube <- o$irf
  if (K > r) {
    cube[, seq_len(K - r), ] <- 0
  }
  if (!is.null(id$size)) {
    for (j in seq_len(r)) {
      cube <- cpp_normalise(
        cube,
        H1,
        id$target[j] - 1L,
        K - r + j - 1L,
        id$size[j]
      )
    }
  }
  list(
    irf = as_irf4(cube, K, H1),
    draw = keep[o$draw],
    accept_rate = sum(o$nacc) / o$tries
  )
}

# ---- (B) Giacomini-Kitagawa-Read robust Bayes --------------------------------

#' Robust-Bayes inference for set-identified proxy-SVARs (Giacomini, Kitagawa & Read)
#'
#' Prior-robust Bayesian inference on impulse responses in a proxy-SVAR identified by
#' exogeneity of `r` proxies plus (optional) sign restrictions, following Giacomini, Kitagawa &
#' Read (2022, *J. Econometrics*) who build on Giacomini & Kitagawa (2021, *Econometrica*).
#' Write the SVAR in orthogonal-reduced form \eqn{y_t=Bx_t+\Sigma_{tr}Q\varepsilon_t} and
#' the first stage \eqn{m_t=Dy_t+Gx_t+\sum_l H_lm_{t-l}+v_t}. The reduced-form parameters
#' \eqn{\phi=(B,\Sigma,D,G,H,\Upsilon)} are updated by the data; the rotation \eqn{Q} is not. The identifying
#' restrictions are the exogeneity restrictions \eqn{D\Sigma_{tr}q_j=0} for the first \eqn{K-r} shocks
#' (the proxies are correlated with the **last `r` shocks**), the user's sign restrictions
#' \eqn{S(\phi,Q)\ge0} and the sign normalisation \eqn{\mathrm{diag}(Q'\Sigma_{tr}^{-1})\ge0}. Instead
#' of one prior for \eqn{Q|\phi}, GKR consider *all* priors supported on the set \eqn{\mathcal{Q}(\phi)}
#' of admissible rotations; this gives, for each \eqn{\phi}, bounds
#' \eqn{[l(\phi),u(\phi)]} of the response, and
#' \itemize{
#'   \item the **set of posterior means** \eqn{[E_\phi l(\phi),\,E_\phi u(\phi)]}, a consistent estimator of the
#'     identified set under strong proxies (GKR sec. 2.4), averaged over draws with a non-empty set;
#'   \item the **robust credible region** of level \eqn{\alpha}: the shortest interval that has posterior
#'     probability \eqn{\ge\alpha} under every prior in the class, centred at
#'     \eqn{\arg\min_\eta \hat z_\alpha(\eta)} with radius \eqn{\min_\eta\hat z_\alpha(\eta)}, where
#'     \eqn{\hat z_\alpha(\eta)} is the \eqn{\alpha}-quantile over draws of
#'     \eqn{\max\{|\eta-l|,|\eta-u|\}};
#'   \item the posterior probability that the identified set is non-empty (plausibility of the sign restrictions).
#' }
#' The bounds are approximated with GKR's **Algorithm 2** (Steps 2.1-2.3 and 3'): per
#' reduced-form draw, up to `nq` rotations satisfying all restrictions are drawn and the
#' extremes taken. This underestimates the identified set at each draw (converging as `nq`
#' grows); GKR's constrained-optimiser Algorithm 1 is **not** implemented. With `size` set, Algorithm 3 is
#' used (responses to a unit shock to `target`; bounds may be unbounded if the impact on `target`
#' can be near zero). `lower`/`upper` of a point-identified shock (e.g. the proxy shock when `r = 1`)
#' coincide, and the set of posterior means then equals the posterior mean of the point-identified response.
#' The single-prior posterior implied by the first accepted rotation of each draw (the conventional
#' Bayesian answer, GKR recommend reporting it too) is returned as `draws`/`bands`.
#'
#' The reduced-form prior is the one in `fit` for \eqn{(B,\Sigma)} (e.g. Minnesota) and Jeffreys
#' for the first stage (GKR sec. 4.1.1 use independent Jeffreys priors; the NIW prior of [bvarm()] also
#' fits their framework). Frequentist coverage is claimed by GKR only for strong proxies
#' (their sec. 3: with a weak proxy the Bernstein-von Mises property fails and the set of posterior means need
#' not be consistent). Check [proxy_fstat()]. The conditioning number of \eqn{D\Sigma_{tr}} diagnostic of
#' their sec. 4.1.6 is not computed.
#'
#' @inheritParams varm_proxy_arw
#' @param signs sign restrictions on all shocks: `K x K` matrix or `K x K x H` array of `+1/-1/0`
#'   (rows variables, columns shocks; the last `r` columns are the proxy shocks), or `NULL`.
#' @param shock column (shock) whose responses are bounded; default `K` (the last proxy shock).
#' @param size `NULL` for responses to a one-standard-deviation shock (GKR default), or the impact on `target`
#'   (single value) for responses to a unit shock (Algorithm 3).
#' @param level credibility levels of the robust credible regions.
#' @param nq accepted rotations per reduced-form draw used for the bounds.
#' @param maxtry rotation draws attempted per reduced-form draw before declaring the set empty.
#' @return a `varm_proxy_gkr` object: `lower`, `upper` (set of posterior means, `(H+1) x K`),
#'   `region` (`(H+1) x K x 2*length(level)`, robust credible region limits ordered like
#'   `bands`), `prob_nonempty`, `lu` (list `l`, `u` of `(H+1) x K x M` arrays over
#'   draws with non-empty set), `draws`/`bands` (single-prior posterior, `shock` column only),
#'   and metadata.
#' @references Giacomini, R., Kitagawa, T. and Read, M. (2022). Robust Bayesian inference in proxy
#'   SVARs. *Journal of Econometrics* 228(1), 107-126.
#'
#'   Giacomini, R. and Kitagawa, T. (2021). Robust Bayesian inference for set-identified models.
#'   *Econometrica* 89(4), 1519-1556.
#' @seealso [varm_proxy_arw()], [varm_proxy_robust()] (frequentist weak-instrument-robust sets)
#' @examples
#' \donttest{
#' set.seed(2)
#' e <- matrix(rnorm(900), 300)
#' y <- matrix(0, 300, 3)
#' for (t in 2:300) y[t, ] <- 0.4 * y[t - 1, ] + e[t, ] %*% t(matrix(c(1, .3, .2, 0, 1, .4, .5, .2, 1), 3))
#' bf <- bvarm(y, p = 1, ndraw = 200, nburn = 100)
#' S <- matrix(0, 3, 3); S[1, 1] <- 1; S[2, 1] <- 1   # shock 1 raises variables 1 and 2 on impact
#' varm_proxy_gkr(bf, m = e[, 3] + rnorm(300, sd = 0.5), signs = S, shock = 1, horizon = 6, nq = 50)
#' }
#' @export
varm_proxy_gkr <- function(
  fit,
  m,
  signs = NULL,
  shock = NULL,
  target = 1L,
  size = NULL,
  pm = 0L,
  horizon = 20L,
  level = c(0.68, 0.90),
  ndraw = NULL,
  nq = 100L,
  maxtry = 2000L,
  sign_horizon = 1L
) {
  proxy_draws_check(fit)
  K <- fit$K
  p <- fit$p
  H <- as.integer(horizon)
  fs <- proxy_first_stage(fit, m, pm)
  r <- fs$r
  if (r >= K) {
    stop("need fewer instruments than variables", call. = FALSE)
  }
  js <- if (is.null(shock)) K else resolve_var(shock, fit$varnames, "shock")
  if (length(js) != 1L || js < 1L || js > K) {
    stop("`shock` must be one column in 1..K", call. = FALSE)
  }
  tg <- -1L
  if (!is.null(size)) {
    tg <- resolve_var(target, fit$varnames, "target")
    if (length(tg) != 1L || length(size) != 1L) {
      stop("with `size`, give one `target` and one `size`", call. = FALSE)
    }
    tg <- tg - 1L
  }
  S <- proxy_sign_cube(signs, K, r, sign_horizon, FALSE)
  Bs <- fit$draws$B
  Ss <- fit$draws$Sigma
  N <- dim(Bs)[3]
  idx <- if (is.null(ndraw) || ndraw >= N) {
    seq_len(N)
  } else {
    sort(sample.int(N, ndraw))
  }
  fd <- proxy_first_stage_draw(fs, length(idx))
  o <- cpp_proxy_rot(
    Bs[,, idx, drop = FALSE],
    Ss[,, idx, drop = FALSE],
    fd$D,
    p,
    H,
    S,
    as.integer(nq),
    1L,
    as.integer(maxtry),
    js - 1L,
    tg,
    if (is.null(size)) 0 else size
  )
  ok <- o$nacc > 0L
  M <- sum(ok)
  if (M == 0L) {
    stop(
      "identified set empty for every draw (sign restrictions incompatible with the data?)",
      call. = FALSE
    )
  }
  l <- o$lower[,, ok, drop = FALSE]
  u <- o$upper[,, ok, drop = FALSE]
  vn <- fit$varnames
  dn <- list(0:H, vn)
  probs <- level_probs(level)
  nl <- length(level)
  alpha <- level
  rk <- match(seq_len(nl), order(level, decreasing = TRUE)) # column of the lower limit
  region <- array(NA_real_, c(H + 1L, K, 2L * nl))
  for (h in seq_len(H + 1L)) {
    for (i in seq_len(K)) {
      for (a in seq_len(nl)) {
        z <- gkr_region(l[h, i, ], u[h, i, ], alpha[a])
        region[h, i, c(rk[a], 2L * nl + 1L - rk[a])] <- c(
          z[1] - z[2],
          z[1] + z[2]
        )
      }
    }
  }
  dimnames(region) <- list(0:H, vn, format(probs))
  # single-prior (first accepted rotation) posterior of the same response
  cube <- o$irf
  if (tg >= 0L) {
    cube <- cpp_normalise(cube, H + 1L, tg, js - 1L, size)
  }
  sp <- as_irf4(cube, K, H + 1L)[, js, , , drop = FALSE]
  sd <- aperm(array(sp, c(K, H + 1L, dim(sp)[4])), c(2, 1, 3))
  out <- list(
    horizon = H,
    p = p,
    K = K,
    r = r,
    varnames = vn,
    shock = js,
    target = if (tg >= 0L) tg + 1L else NA_integer_,
    size = size,
    level = level,
    method = "gkr",
    nq = nq,
    ndraw = length(idx),
    lower = structure(apply(l, 1:2, mean), dimnames = dn),
    upper = structure(apply(u, 1:2, mean), dimnames = dn),
    region = region,
    prob_nonempty = M / length(idx),
    nacc = o$nacc,
    lu = list(l = l, u = u),
    estimate = structure(apply(sd, 1:2, stats::median), dimnames = dn)
  )
  out <- proxy_bands(
    out,
    array(sd, c(H + 1L, K, 1L, dim(sd)[3])),
    level,
    vn,
    NULL
  )
  structure(out, class = "varm_proxy_gkr")
}

# smallest robust credible interval for one scalar: centre and radius (GKR Algorithm 1, Step 5)
gkr_region <- function(l, u, alpha) {
  f <- function(eta) {
    stats::quantile(pmax(abs(eta - l), abs(eta - u)), alpha, names = FALSE)
  }
  if (all(l == u) && length(unique(l)) == 1L) {
    return(c(l[1], 0))
  }
  lo <- min(l)
  hi <- max(u)
  if (lo == hi) {
    return(c(lo, 0))
  }
  g <- seq(lo, hi, length.out = 101L)
  v <- vapply(g, f, 0)
  k <- which.min(v)
  a <- g[max(1L, k - 1L)]
  b <- g[min(101L, k + 1L)]
  if (!(a < b)) {
    return(c(g[k], v[k]))
  }
  op <- stats::optimize(f, c(a, b))
  if (op$objective < v[k]) c(op$minimum, op$objective) else c(g[k], v[k])
}

#' @export
print.varm_proxy_gkr <- function(x, hs = c(0L, 4L, 8L, 12L), digits = 3, ...) {
  cat(
    "<varm_proxy_gkr> robust-Bayes bounds, shock ",
    x$shock,
    " (",
    x$r,
    " proxy",
    if (x$r > 1) "ies",
    "; ",
    if (is.null(x$size)) {
      "1 s.d. shock"
    } else {
      paste0("unit shock, size ", x$size)
    },
    "), horizon 0..",
    x$horizon,
    "\n",
    x$ndraw,
    " reduced-form draws, nq = ",
    x$nq,
    "; P(identified set non-empty) = ",
    format(x$prob_nonempty, digits = 3),
    "\n",
    sep = ""
  )
  hs <- hs[hs <= x$horizon]
  nl <- length(x$level)
  cat(
    "set of posterior means [lower, upper] and ",
    format(100 * x$level[nl]),
    "% robust credible region:\n",
    sep = ""
  )
  for (i in seq_len(x$K)) {
    cat(x$varnames[i], ":\n", sep = "")
    tab <- cbind(
      lower = x$lower[hs + 1L, i],
      upper = x$upper[hs + 1L, i],
      rcr_lo = x$region[hs + 1L, i, 1L],
      rcr_hi = x$region[hs + 1L, i, 2L * nl]
    )
    rownames(tab) <- paste0("h=", hs)
    print(round(tab, digits))
  }
  invisible(x)
}
