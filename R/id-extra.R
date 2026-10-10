# Extra identification schemes: proxy-SVAR, zero+sign restrictions, A-B restrictions, heteroskedasticity.
# Each id_* object here carries a `structural(Bs, Ss, K, p, H, id, varnames)` closure that
# structural_irf() calls (one dispatch line in R/id.R).

#' External-instrument (proxy) SVAR identification
#'
#' \code{id_proxy()} identifies a single structural shock \eqn{\varepsilon_{1t}} from an external
#' instrument \eqn{z_t} (Mertens & Ravn 2013; Stock & Watson 2018; Gertler & Karadi 2015).
#' With \eqn{u_t = s\,\varepsilon_{1t} + \dots}, the instrument satisfies relevance
#' \eqn{E[z_t\varepsilon_{1t}]=\phi\neq0} and exogeneity \eqn{E[z_t\varepsilon_{jt}]=0,\ j>1}, hence
#' \deqn{\gamma = E[u_t z_t] = s\,\phi,\qquad s \propto \gamma .}
#' The impact vector is identified up to scale; it is normalised so that the impact response of
#' variable `target` equals `size`: \eqn{s = \text{size}\cdot\gamma/\gamma_{target}}. The response
#' at horizon \eqn{h} is \eqn{\Phi_h s}, \eqn{\Phi_h} the reduced-form MA coefficients. The estimate is
#' the first stage of the residual on the proxy: \eqn{\hat\gamma = T^{-1}\sum_t \hat u_t z_t}, and
#' relevance is judged by the first-stage F of the `target` residual on \eqn{z} (see [proxy_fstat()]).
#'
#' Only column `shock` of the structural impact matrix is identified: the other columns of the
#' returned IRF array are zero, so FEVD of the other shocks is meaningless (the FEVD column of the
#' proxy shock is correct only if `size` is the shock standard deviation, which a proxy cannot deliver).
#'
#' Uncertainty needs the data, not only \eqn{(B,\Sigma)}: use [proxy_draws()] which bootstraps
#' residuals and proxy jointly and carries one \eqn{\gamma} per draw. Using `irf(fit, id_proxy(z, fit))`
#' with `ci != "none"` holds \eqn{\gamma} fixed and therefore understates uncertainty.
#'
#' *Economic intuition.* A macroeconomist uses a proxy when a credible narrative/high-frequency
#' measure of one shock is available (monetary-policy surprises around FOMC announcements, tax
#' narratives, oil supply news) but a full recursive/sign scheme is not trusted. Pitfalls: weak
#' instruments (F below about 10 makes the normalised IRF very wide and non-normal; use the
#' Montiel Olea-Stock-Watson F_hac or Anderson-Rubin bands, see [varm_proxy_robust()]), proxy
#' contaminated by other shocks (exogeneity), and the sample used for the proxy being shorter than the VAR
#' (set missing proxy values to 0 only if the instrument is truly 0 there).
#'
#' @param z proxy: numeric vector of length `nrow(fit$Y)` (aligned with the data; first `p` values dropped).
#' @param fit a [varm()] fit (used to compute the point `gamma`).
#' @param gamma alternatively, a `K` vector or `K x N` matrix of reduced-form impact covariances (one
#'   column per reduced-form draw, as produced by [proxy_draws()]).
#' @param target <[`tidy-select`][tidyselect::language]> variable whose impact response is normalised
#'   (name or position).
#' @param size impact response of `target` to the shock (e.g. 0.25 for a 25bp surprise).
#' @param shock <[`tidy-select`][tidyselect::language]> column of the IRF array that holds the proxy shock
#'   (default `target`).
#' @return an object of class `varm_id`.
#' @seealso [proxy_draws()], [varm_proxy()]
#' @examples
#' fit <- varm(yraw, p = 2)
#' z <- c(rep(0, 2), resid(fit)[, "rate"] + rnorm(nrow(resid(fit))))
#' id <- id_proxy(z, fit, target = rate)
#' irf(fit, id, horizon = 4, ci = "none")
#' @family identification
#' @export
id_proxy <- function(
  z = NULL,
  fit = NULL,
  gamma = NULL,
  target = 1L,
  size = 1,
  shock = target
) {
  target <- rlang::enquo(target)
  shock <- if (missing(shock)) target else rlang::enquo(shock)
  fstat <- NULL
  if (is.null(gamma)) {
    if (is.null(z) || is.null(fit)) {
      rlang::abort("supply `gamma`, or both `z` and `fit`")
    }
    z <- check_proxy(z, fit)
    gamma <- drop(crossprod(fit$U, z[-seq_len(fit$p)])) / nrow(fit$U)
    fstat <- proxy_fstat(
      fit,
      z,
      select_one(target, fit$varnames, "target")
    )
  }
  new_id(
    "proxy",
    gamma = as.matrix(gamma),
    target = target,
    size = size,
    shock = shock,
    fstat = fstat,
    structural = proxy_structural
  )
}

check_proxy <- function(z, fit) {
  z <- as.numeric(z)
  if (length(z) != nrow(fit$Y) || anyNA(z)) {
    rlang::abort("`z` must have no NAs and length nrow(fit$Y)")
  }
  z
}

proxy_structural <- function(Bs, Ss, K, p, H, id, varnames) {
  G <- id$gamma
  if (nrow(G) != K) {
    rlang::abort("`gamma` must have K rows")
  }
  if (ncol(G) != 1L && ncol(G) != dim(Bs)[3]) {
    rlang::abort("`gamma` must have 1 or N columns")
  }
  j <- select_one(id$target, varnames, "target")
  s <- select_one(id$shock, varnames, "shock")
  G <- id$size * sweep(G, 2, G[j, ], "/")
  cube <- cpp_proxy_irf(Bs, G, p, H, s - 1L)
  list(
    irf = as_irf4(cube, K, H + 1L),
    draw = seq_len(dim(Bs)[3]),
    accept_rate = NA_real_
  )
}

#' Joint residual-proxy bootstrap for proxy SVARs
#'
#' Draws `nboot` bootstrap samples of \eqn{(B,\Sigma,\gamma)}: `"wild"` multiplies \eqn{(\hat u_t, z_t)}
#' by a common Rademacher weight (Mertens-Ravn style; NOT recommended, since \eqn{w_t^2=1} leaves
#' \eqn{T^{-1}\sum u_t z_t} unchanged, so \eqn{\gamma} only varies through re-estimated VAR
#' coefficients and bands are far too narrow, as shown by Jentsch & Lunsford 2019 and in
#' `validation/ident-extra.R`); `"mbb"` (default) resamples moving blocks of the pairs \eqn{(\hat u_t, z_t)}
#' (block length `block`, default the Jentsch-Lunsford rule \eqn{\mathrm{round}(5.03\,T^{1/4})}), preserving serial dependence in
#' both (Jentsch & Lunsford 2019). Each draw re-estimates the VAR on the bootstrap data and
#' recomputes \eqn{\gamma^* = T^{-1}\sum u^*_t z^*_t}. Percentile bands of the normalised IRF
#' are only reliable with a strong instrument.
#' @inheritParams id_proxy
#' @param nboot number of draws.
#' @param type `"mbb"` (default, valid) or `"wild"` (invalid for gamma uncertainty, kept for comparison).
#' @param block block length for `"mbb"`.
#' @return a list with `draws` (a [draws] object, plus `gamma`), `id` (an [id_proxy()] carrying one
#'   gamma per draw) and `fstat`. Use `irf(res$draws, res$id)`.
#' @examples
#' fit <- varm(yraw, p = 2)
#' z <- c(rep(0, 2), resid(fit)[, "rate"] + rnorm(nrow(resid(fit))))
#' pd <- proxy_draws(fit, z, nboot = 20, target = rate)
#' irf(pd$draws, pd$id, horizon = 4)
#' @family identification
#' @export
proxy_draws <- function(
  fit,
  z,
  nboot = 1000L,
  type = c("mbb", "wild"),
  block = NULL,
  target = 1L,
  size = 1,
  shock = target
) {
  type <- match.arg(type)
  target <- rlang::enquo(target)
  shock <- if (missing(shock)) target else rlang::enquo(shock)
  z <- check_proxy(z, fit)
  block <- block %||% round(5.03 * nrow(fit$U)^(1 / 4))
  b <- cpp_proxy_boot(
    fit$Y,
    fit$p,
    fit$det,
    fit$B,
    fit$U,
    z[-seq_len(fit$p)],
    as.integer(nboot),
    c(wild = 0L, mbb = 1L)[[type]],
    as.integer(block)
  )
  d <- structure(
    list(
      B = b$B,
      Sigma = b$Sigma,
      p = fit$p,
      K = fit$K,
      varnames = fit$varnames,
      method = paste0("proxy_", type),
      fit = fit,
      gamma = b$gamma
    ),
    class = "varm_draws"
  )
  list(
    draws = d,
    id = id_proxy(
      gamma = b$gamma,
      target = !!target,
      size = size,
      shock = !!shock
    ),
    fstat = proxy_fstat(fit, z, select_one(target, fit$varnames, "target"))
  )
}

zs_array <- function(x, K, horizon, what) {
  if (is.null(x)) {
    return(array(0, c(K, K, 1L)))
  }
  if (length(dim(x)) == 2L) {
    x <- array(rep(x, horizon), c(dim(x), horizon))
  }
  if (length(dim(x)) != 3L || any(dim(x)[1:2] != K)) {
    rlang::abort(paste0("`", what, "` must be K x K or K x K x H"))
  }
  x[is.na(x)] <- 0
  x
}

#' Sign restrictions combined with zero restrictions (Arias, Rubio-Ramirez & Waggoner 2018)
#'
#' Let \eqn{\Theta_h = \Phi_h L Q} be the structural IRF at horizon \eqn{h}, with \eqn{LL'=\Sigma} and
#' \eqn{Q} orthogonal, and \eqn{\Theta_\infty = (I-\sum_l A_l)^{-1} L Q}. A zero restriction on the response
#' of variable \eqn{i} to shock \eqn{j} at horizon \eqn{h} reads \eqn{e_i'\Phi_h L\, q_j = 0}, linear in the
#' column \eqn{q_j}. Following ARRW Algorithm 2, shocks are processed in decreasing order of
#' their number of zero restrictions \eqn{z_j}; column \eqn{q_j} is drawn uniformly on the unit sphere of
#' the null space of \eqn{[Z_j L;\ q_1';\dots;q_{j-1}']}, which exists iff \eqn{z_{(j)} \le K-j}. The
#' draw is kept if all sign restrictions (up to a per-shock sign flip) hold. Zeros hold exactly by
#' construction and \eqn{PP'=\Sigma}.
#'
#' **Missing importance weights.** ARRW additionally reweight draws by a Jacobian term
#' (\eqn{|\det A_0|^{-(2n+m+1)}} over \eqn{\prod_i |\det V_i|}) so that they are uniform on the
#' restricted set under the conjugate prior. These weights are NOT implemented here: accepted draws are
#' *unweighted* samples from the null-space sampler, which is the exact ARRW posterior only when no
#' zero restrictions are present (then it coincides with [id_sign()]) or when the zero restrictions
#' fully pin down the rotation up to sign (then the weights are constant across rotations, though they
#' still vary over the reduced-form draws). Treat bands under partial zero-restriction identification
#' as approximate.
#'
#' *Economic intuition.* Use it when theory gives a few exact zeros (e.g. monetary policy has no
#' contemporaneous effect on output and prices, or a demand shock has no long-run effect on output)
#' and signs for the rest. Pitfalls: bands mix identification uncertainty and sampling
#' uncertainty; the implied prior over impulse responses is informative; zeros on many shocks can
#' make the problem exactly identified (sign information then irrelevant).
#'
#' @param sign `K x K` matrix or `K x K x H` array of `+1/-1/0` as in [id_sign()] (`NULL` for none).
#' @param zero `K x K` matrix or `K x K x H` array with `1` where the response (variable x shock) is
#'   zero at that horizon (impact only when a matrix: `zero_horizon = 1`).
#' @param zero_lr `K x K` 0/1 matrix, `1` = long-run (cumulative) response is zero.
#' @param horizon horizon count over which a *matrix* `sign` is imposed.
#' @param zero_horizon horizon count for a *matrix* `zero` (default 1: impact only).
#' @param nrot,maxtry as in [id_sign()].
#' @return a `varm_id` object (usable in [irf()] on `varm_draws`/`varm`).
#' @seealso [id_sign()], [id_restrictions()], [irf_fry_pagan()]
#' @examples
#' S <- matrix(0, 3, 3)
#' S[1, 1] <- 1
#' S[3, 3] <- 1
#' Z <- matrix(0, 3, 3)
#' Z[3, 2] <- 1
#' id <- id_zerosign(S, Z)
#' irf(draws(varm(yraw, p = 2), "mc", 20), id, horizon = 4)
#' @family identification
#' @export
id_zerosign <- function(
  sign = NULL,
  zero = NULL,
  zero_lr = NULL,
  horizon = 1L,
  zero_horizon = 1L,
  nrot = 1L,
  maxtry = 1000L
) {
  K <- dim(sign %||% zero %||% zero_lr)[1]
  S <- zs_array(sign, K, horizon, "sign")
  Z <- zs_array(zero, K, zero_horizon, "zero")
  if (!all(S %in% c(-1, 0, 1))) {
    rlang::abort("sign restrictions must be -1, 0, +1")
  }
  lr <- zero_lr %||% matrix(0, K, K)
  lr[is.na(lr)] <- 0
  if (!all(Z %in% 0:1) || !all(lr %in% 0:1)) {
    rlang::abort("zero restrictions must be 0/1")
  }
  new_id(
    "zerosign",
    S = S,
    Z = Z,
    lr = lr,
    nrot = as.integer(nrot),
    maxtry = as.integer(maxtry),
    structural = zs_structural
  )
}

zs_structural <- function(Bs, Ss, K, p, H, id, varnames) {
  if (dim(id$S)[1] != K) {
    rlang::abort(
      "restrictions dimension does not match the number of variables"
    )
  }
  out <- cpp_zs_irf(Bs, Ss, p, H, id$S, id$Z, id$lr, id$nrot, id$maxtry)
  if (out$accepted == 0L) {
    rlang::abort(
      paste0("no draw satisfied the restrictions (tries = ", out$tries, ")")
    )
  }
  list(
    irf = as_irf4(out$irf, K, H + 1L),
    draw = out$draw,
    P = out$P,
    accept_rate = out$accepted / out$tries
  )
}

#' Zero (short-run and long-run) restrictions in an A-B / Blanchard-Quah style
#'
#' General exactly- or under-identifying zero pattern on the impact matrix \eqn{P} (\eqn{u_t=P\varepsilon_t})
#' and on the long-run matrix \eqn{\Theta_\infty = (I-\sum A_l)^{-1}P}. Entries equal to `0` are
#' restricted to zero, `NA` (or anything else) are free. For a just-identified pattern (after ordering
#' shocks by number of zeros the counts are \eqn{K-1,K-2,\dots,0}) the rotation is unique up to
#' column signs, found in closed form by sequential null-space projections
#' (the ARRW algorithm of [id_zerosign()], no numerical optimiser needed), so \eqn{\Sigma=PP'}
#' and every zero hold to machine precision. Shock signs are normalised so that the impact response
#' of variable `j` to shock `j` (or the first free entry of column `j`) is positive. Over-identified patterns are
#' rejected (more zeros than the counting allows); under-identified patterns give random draws from the
#' admissible set (see the missing-weights note in [id_zerosign()]).
#'
#' *Economic intuition.* The usual structural-VAR workhorse: Sims' recursive scheme is
#' `impact = lower triangular`, Blanchard-Quah is `longrun = lower triangular`, and mixed schemes
#' (e.g. Gali 1992, Blanchard-Perotti with a long-run block) combine both. Pitfall: the zeros are
#' identifying assumptions that cannot be tested in the just-identified case.
#' @param impact,longrun `K x K` matrices with `0` for zero restrictions and `NA` for free elements.
#' @return a `varm_id` object.
#' @examples
#' lower <- matrix(NA, 3, 3)
#' lower[upper.tri(lower)] <- 0
#' irf(varm(yraw, p = 2), id_restrictions(impact = lower), 3, ci = "none")
#' @family identification
#' @export
id_restrictions <- function(impact = NULL, longrun = NULL) {
  K <- nrow(impact %||% longrun)
  Z <- if (!is.null(impact)) {
    array(as.numeric(!is.na(impact) & impact == 0), c(K, K, 1L))
  }
  lr <- if (!is.null(longrun)) 1 * (!is.na(longrun) & longrun == 0)
  # one positive impact per shock fixes its orientation: the own-variable
  # response when it is unrestricted, else the first unrestricted response
  oriented <- function(j) {
    if (is.null(impact)) {
      return(j)
    }
    free <- which(is.na(impact[, j]) | impact[, j] != 0)
    if (j %in% free) j else free[1]
  }
  S <- array(0, c(K, K, 1L))
  for (j in seq_len(K)) {
    S[oriented(j), j, 1L] <- 1
  }
  id <- id_zerosign(sign = S, zero = Z, zero_lr = lr, maxtry = 10L)
  id$type <- "restrictions"
  class(id) <- c("id_restrictions", "varm_id")
  id
}

#' Identification through heteroskedasticity (Rigobon 2003), two regimes
#'
#' With two covariance regimes \eqn{\Sigma_1=PP'} and \eqn{\Sigma_2=P\Lambda P'},
#' \eqn{\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_K)} with distinct elements, \eqn{P} is identified
#' up to column order and sign: with \eqn{\Sigma_1=L_1L_1'}, take the eigendecomposition
#' \eqn{L_1^{-1}\Sigma_2L_1^{-T}=V\Lambda V'}, \eqn{P=L_1V}. Shocks are ordered by decreasing
#' \eqn{\lambda_j} (relative variance in regime 2), signs chosen so \eqn{\mathrm{diag}(P)>0}.
#' The result is a *point* impact matrix (\eqn{\Sigma} from the full sample need not equal
#' \eqn{PP'}; \eqn{P} satisfies it for the regime-1 covariance), returned as an [id_given()] object with
#' attribute `lambda`. It does not propagate sampling uncertainty (all draws share one P).
#'
#' *Economic intuition.* Useful when volatility shifts (Great Moderation, crisis episodes,
#' announcement days) are plausibly shock-specific, with no need for any economic restriction;
#' the shocks are statistical and must be labelled afterwards. Pitfalls: needs
#' distinct variance ratios (weak identification otherwise), regimes known/exogenous, and
#' constant \eqn{P} across regimes.
#' @param fit a [varm()] fit.
#' @param regime logical or 0/1/2 vector of length `nrow(fit$Y)`, `TRUE`/`2` = regime 2.
#' @examples
#' fit <- varm(yraw, p = 1)
#' regime <- rep(0:1, length.out = nrow(yraw))
#' id_hetero(fit, regime)
#' @family identification
#' @export
id_hetero <- function(fit, regime) {
  r <- as.numeric(regime)[-seq_len(fit$p)]
  if (length(r) != nrow(fit$U)) {
    rlang::abort("`regime` must have length nrow(fit$Y)")
  }
  hi <- r == max(r)
  cov_regime <- function(sel) crossprod(fit$U[sel, , drop = FALSE]) / sum(sel)
  L1 <- t(chol(cov_regime(!hi)))
  e <- eigen(solve(L1, t(solve(L1, cov_regime(hi)))), symmetric = TRUE)
  P <- L1 %*% e$vectors
  P <- sweep(P, 2, ifelse(diag(P) < 0, -1, 1), "*")
  dimnames(P) <- list(fit$varnames, fit$varnames)
  id <- id_given(P)
  attr(id, "lambda") <- e$values
  id
}

#' Fry-Pagan (2011) median-target rotation
#'
#' Among the accepted draws of a set-identified [irf()] (sign / zero-sign), picks the single draw
#' whose IRFs are closest, in standardised squared distance over all responses, shocks and horizons,
#' to the pointwise median. It returns one internally consistent structural model instead of the
#' incoherent pointwise median.
#'
#' *Economic intuition.* Report this alongside the set bands when you need one coherent
#' shock series, historical decomposition or FEVD. Pitfall: the median-target is itself a function
#' of the prior over rotations and is not "the" estimate.
#' @param x a `varm_irf` object with draws (from sign-type identification).
#' @return list with `index` (into `x$draws[,,,index]` and `x$draw`), `irf` (`[K,K,H+1]`),
#'   `impact` (`K x K` impact matrix) and `varm_draw` (index of the reduced-form draw it came from).
#' @examples
#' S <- matrix(c(1, 1, 1, 1, -1, 0, 1, 0, 0), 3, 3, byrow = TRUE)
#' r <- irf(draws(varm(yraw, p = 2), "mc", 30), id_sign(S), horizon = 4)
#' irf_fry_pagan(r)$index
#' @family identification
#' @export
irf_fry_pagan <- function(x) {
  i <- fry_pagan(x$draws)
  a <- x$draws[,,, i]
  dimnames(a) <- dimnames(x$estimate)
  list(index = i, irf = a, impact = a[,, 1], varm_draw = x$draw[i])
}
