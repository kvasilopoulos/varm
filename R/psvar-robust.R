# Internal: lagged regressor matrix (T-p) x m, same column order as fit$B.
proxy_X <- function(fit) {
  p <- fit$p
  Tn <- nrow(fit$Y)
  cbind(
    do.call(
      cbind,
      lapply(seq_len(p), function(l) {
        fit$Y[(p + 1 - l):(Tn - l), , drop = FALSE]
      })
    ),
    fit$det[-seq_len(p), , drop = FALSE]
  )
}

# Internal: Gamma = T^-1 sum u_t m_t and the HAC asymptotic covariance W of
# sqrt(T) * (vec(A), Gamma), A = [A_1 ... A_p] (K x Kp), via the stacked influence
# functions of the OLS coefficients and of Gamma-hat (Newey-West, Bartlett weights).
proxy_cov <- function(fit, me, nwlags = NULL) {
  K <- fit$K
  p <- fit$p
  U <- fit$U
  Tp <- nrow(U)
  if (is.null(nwlags)) {
    nwlags <- floor(4 * (Tp / 100)^(2 / 9))
  }
  X <- proxy_X(fit)
  Q <- crossprod(X) / Tp
  Xi <- X %*% solve(Q) # row t: (Q^-1 x_t)'
  g <- drop(crossprod(U, me)) / Tp
  qz <- drop(crossprod(X, me)) / Tp # E[x z]
  # psi_B[t, , ] = Q^-1 x_t u_t' (m x K); Gamma influence: u z - Gamma - psi_B' qz
  psiA <- matrix(0, Tp, K * K * p)
  psiG <- matrix(0, Tp, K)
  for (t in seq_len(Tp)) {
    PB <- tcrossprod(Xi[t, ], U[t, ])
    psiA[t, ] <- as.vector(t(PB[seq_len(K * p), , drop = FALSE]))
    psiG[t, ] <- U[t, ] * me[t] - g - drop(crossprod(PB, qz))
  }
  S <- scale(cbind(psiA, psiG), center = TRUE, scale = FALSE)
  W <- crossprod(S) / Tp
  for (l in seq_len(nwlags)) {
    Gl <- crossprod(
      S[-seq_len(l), , drop = FALSE],
      S[seq_len(Tp - l), , drop = FALSE]
    ) /
      Tp
    W <- W + (1 - l / (nwlags + 1)) * (Gl + t(Gl))
  }
  list(
    Gamma = g,
    W = W,
    iA = seq_len(K * K * p),
    iG = K * K * p + seq_len(K),
    Tp = Tp,
    nwlags = nwlags
  )
}

# Internal: d vec(Phi_h) / d vec(A) for h = 0..H (list of K^2 x K^2 p matrices), Luetkepohl (2005, Prop. 3.6).
proxy_dphi <- function(B, K, p, H) {
  Phi <- cpp_ma(B, K, p, H)
  A <- t(B[seq_len(K * p), , drop = FALSE])
  Ac <- rbind(A, cbind(diag(K * (p - 1)), matrix(0, K * (p - 1), K)))[
    seq_len(K * p),
    ,
    drop = FALSE
  ]
  if (p == 1L) {
    Ac <- A
  }
  J <- cbind(diag(K), matrix(0, K, K * (p - 1)))
  pw <- vector("list", H + 1L)
  pw[[1]] <- diag(K * p)
  for (i in seq_len(H)) {
    pw[[i + 1]] <- pw[[i]] %*% t(Ac)
  }
  G <- lapply(0:H, function(h) {
    out <- matrix(0, K * K, K * K * p)
    for (m in seq_len(h) - 1L) {
      if (h > 0) out <- out + kronecker(J %*% pw[[h - m]], Phi[,, m + 1L])
    }
    out
  })
  list(Phi = Phi, G = G)
}

#' Weak-instrument-robust confidence sets for proxy-SVAR impulse responses
#'
#' Montiel Olea, Stock & Watson (2021) Anderson-Rubin-type confidence sets for
#' the impulse responses to a shock identified with a single external instrument.
#' For each variable \eqn{i} and horizon \eqn{h}, the set collects all \eqn{\lambda}
#' for which the (delta-method, HAC) Wald test of
#' \eqn{s\, e_i'\Phi_h\Gamma - \lambda\Gamma_j = 0} does not reject, where
#' \eqn{\Gamma = E[u_t m_t]}, \eqn{j} is the `target` variable and \eqn{s} is `size`.
#' The test is a quadratic \eqn{a\lambda^2 + b\lambda + c \le 0}, so the set is one of
#' four types (reported in `case`):
#' \describe{
#'   \item{1}{`a > 0`, discriminant `> 0`: bounded interval `[lower, upper]`.}
#'   \item{2}{`a < 0`, discriminant `> 0`: the *complement* of an interval, i.e.
#'     `(-Inf, r1] U [r2, Inf)`. `lower = -Inf`, `upper = Inf` and the excluded
#'     interval `(r1, r2)` is in `gap`.}
#'   \item{3}{`a > 0`, discriminant `< 0`: empty set (`lower = upper = NA`);
#'     signals misspecification.}
#'   \item{4}{otherwise (`a < 0`, discriminant `<= 0`): the whole real line
#'     (`-Inf`, `Inf`); the instrument is uninformative.}
#'   \item{0}{impact response of `target`, normalised to `size`.}
#' }
#' Unlike the bootstrap bands of [varm_proxy()], these have correct coverage
#' regardless of instrument strength. The covariance of
#' \eqn{(\mathrm{vec}(A), \Gamma)} uses the OLS influence functions with a Bartlett
#' (Newey-West) HAC estimator; derivatives of the MA matrices follow Luetkepohl (2005).
#' @inheritParams varm_proxy
#' @param level confidence levels (in `(0, 1)`).
#' @param nwlags Newey-West truncation lag; default `floor(4 * (T/100)^(2/9))`.
#' @return a `varm_proxy_robust` object: `estimate` (`(H+1) x K` point IRF),
#'   `lower`, `upper`, `case` (arrays `(H+1) x K x length(level)`), `gap`
#'   (`(H+1) x K x length(level) x 2`, excluded interval in case 2), `fstat`
#'   (see [proxy_fstat()]), `level`.
#' @references Montiel Olea, J. L., Stock, J. H. and Watson, M. W. (2021). Inference
#'   in structural vector autoregressions identified with an external instrument.
#'   *Journal of Econometrics* 225(1), 74-87.
#' @examples
#' set.seed(1)
#' e <- matrix(rnorm(600), 300)
#' y <- matrix(0, 300, 2)
#' for (t in 2:300) y[t, ] <- 0.5 * y[t - 1, ] + c(e[t, 1] + 0.3 * e[t, 2], e[t, 2] + 0.5 * e[t, 1])
#' fit <- varm(y, p = 1)
#' r <- varm_proxy_robust(fit, m = e[, 1] + rnorm(300, sd = 0.5), horizon = 6)
#' r$lower[, , 2]
#' @export
varm_proxy_robust <- function(
  fit,
  m,
  target = 1L,
  size = 1,
  horizon = 20L,
  level = c(0.68, 0.90),
  nwlags = NULL
) {
  if (!inherits(fit, "varm")) {
    stop("`fit` must be a varm object", call. = FALSE)
  }
  K <- fit$K
  p <- fit$p
  Tn <- nrow(fit$Y)
  H <- as.integer(horizon)
  m <- as.numeric(m)
  if (length(m) != Tn || anyNA(m)) {
    stop("`m` must have no NAs and length nrow(fit$Y)", call. = FALSE)
  }
  if (any(level <= 0 | level >= 1)) {
    stop("`level` must be in (0, 1)", call. = FALSE)
  }
  j <- resolve_var(target, fit$varnames, "target")
  cv <- proxy_cov(fit, m[-seq_len(p)], nwlags)
  g <- cv$Gamma
  Tp <- cv$Tp
  WA <- cv$W[cv$iA, cv$iA]
  WAG <- cv$W[cv$iA, cv$iG]
  WG <- cv$W[cv$iG, cv$iG]
  dp <- proxy_dphi(fit$B, K, p, H)
  crit <- stats::qchisq(level, 1)
  L <- length(level)
  dn <- list(0:H, fit$varnames, format(level))
  lower <- upper <- cs <- array(NA_real_, c(H + 1L, K, L), dimnames = dn)
  gap <- array(NA_real_, c(H + 1L, K, L, 2L))
  est <- matrix(0, H + 1L, K, dimnames = dn[1:2])
  ej <- diag(K)[, j]
  for (h in 0:H) {
    Phi <- dp$Phi[,, h + 1L]
    G <- dp$G[[h + 1L]]
    for (i in seq_len(K)) {
      ei <- diag(K)[, i]
      est[h + 1L, i] <- size * sum(Phi[i, ] * g) / g[j]
      f0 <- size * sum(Phi[i, ] * g)
      w <- size * drop(as.vector(outer(ei, g)) %*% G)
      u <- size * Phi[i, ]
      wWv <- sum(w * (WAG %*% ej)) + sum(u * (WG %*% ej))
      cc0 <- drop(w %*% WA %*% w + 2 * w %*% WAG %*% u + u %*% WG %*% u)
      for (l in seq_len(L)) {
        a <- Tp * g[j]^2 - crit[l] * WG[j, j]
        b <- -2 * Tp * f0 * g[j] + 2 * crit[l] * wWv
        c0 <- Tp * f0^2 - crit[l] * cc0
        D <- b^2 - 4 * a * c0
        r <- if (D > 0) {
          sort((-b + c(-1, 1) * sqrt(D)) / (2 * a))
        } else {
          c(NA_real_, NA_real_)
        }
        if (a > 0 && D > 0) {
          cs[h + 1L, i, l] <- 1
          lower[h + 1L, i, l] <- r[1]
          upper[h + 1L, i, l] <- r[2]
        } else if (a < 0 && D > 0) {
          cs[h + 1L, i, l] <- 2
          lower[h + 1L, i, l] <- -Inf
          upper[h + 1L, i, l] <- Inf
          gap[h + 1L, i, l, ] <- r
        } else if (a > 0) {
          cs[h + 1L, i, l] <- 3
        } else {
          cs[h + 1L, i, l] <- 4
          lower[h + 1L, i, l] <- -Inf
          upper[h + 1L, i, l] <- Inf
        }
      }
    }
  }
  # impact response of the target is the normalisation
  lower[1, j, ] <- upper[1, j, ] <- size
  cs[1, j, ] <- 0
  gap[1, j, , ] <- NA
  structure(
    list(
      estimate = est,
      lower = lower,
      upper = upper,
      case = cs,
      gap = gap,
      level = level,
      horizon = H,
      varnames = fit$varnames,
      target = j,
      size = size,
      fstat = proxy_fstat(fit, m, j, nwlags)
    ),
    class = "varm_proxy_robust"
  )
}

#' @export
print.varm_proxy_robust <- function(x, ...) {
  cat(
    "<varm_proxy_robust> Montiel Olea-Stock-Watson sets, shock to ",
    x$varnames[x$target],
    " (size ",
    x$size,
    "), horizon 0..",
    x$horizon,
    "; levels ",
    paste(x$level, collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "HAC first-stage F = ",
    format(x$fstat[["F_hac"]], digits = 4),
    "\n",
    sep = ""
  )
  tab <- table(factor(
    x$case,
    0:4,
    c("normalised", "bounded", "complement", "empty", "real line")
  ))
  print(tab[tab > 0])
  invisible(x)
}
