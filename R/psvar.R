#' Proxy-SVAR (external instrument) identification and bootstrap inference
#'
#' Mertens & Ravn (2013) / Stock & Watson (2012) proxy identification of `r` structural
#' shocks in a reduced-form [varm()] fit, using `r` external instruments. With one instrument
#' the impact vector is \eqn{\gamma = T^{-1}\sum_t u_t m_t} (covariance of the reduced-form
#' residuals with the proxy), scaled so that the impact response of `target` equals `size`.
#' Responses are \eqn{\Theta_h \gamma}, reusing the MA machinery of [irf()]. With `r > 1`
#' instruments the Mertens-Ravn (2013) `B11`/`B21` formulae are used (see Details).
#'
#' @details
#' **Inference (`method`).**
#' \describe{
#'   \item{`"wild"`}{Mertens & Ravn (2013) / Gertler & Karadi (2015)-style recursive-design
#'     wild bootstrap: Rademacher weights multiply the residuals and the proxy jointly.}
#'   \item{`"wild_normal"`}{same with standard-normal weights (the other variant studied
#'     by Jentsch & Lunsford).}
#'   \item{`"mbb"`}{residual-based moving block bootstrap of Jentsch & Lunsford: blocks
#'     of length `block` are drawn jointly from the residuals and the proxy, each is
#'     recentred by the mean of the corresponding within-block position over all
#'     blocks, the VAR is re-simulated recursively and re-estimated. Block length default is
#'     the Jentsch-Lunsford rule \eqn{\ell = 5.03\,T^{1/4}} (\eqn{\ell = 20} at
#'     \eqn{T = 250}).}
#' }
#' Jentsch & Lunsford (2019, 2022) show that the wild bootstrap is asymptotically *invalid* for
#' proxy-SVAR impulse responses, because multiplying residuals and proxy by the same weight
#' leaves the sample covariance \eqn{T^{-1}\sum m_t u_t'} (almost) fixed, so the sampling
#' variability of the impact vector is not reproduced, and that the moving block bootstrap
#' is valid. Mertens & Ravn (2019) dispute the practical conclusions. The wild bootstrap
#' remains available (`method = "wild"`) to reproduce Mertens-Ravn / Gertler-Karadi bands.
#' **The default changed from `"wild"` to `"mbb"`** on the strength of that theory and of the
#' Monte Carlo in `validation/proxy-coverage.R` (wild-bootstrap 90% coverage of the impact
#' response of a non-target variable was about 14% vs 87% for the block bootstrap). Neither bootstrap is robust to weak instruments
#' ([varm_proxy_robust()] is). Intervals are equal-tailed percentiles, as in the papers.
#'
#' **Several instruments.** `m` may be a `T x r` matrix and `target` a vector of `r`
#' variables: `r` instruments identify `r` shocks. As in the Mertens-Ravn algorithm
#' (and the SVAR-IV MATLAB code of Miranda-Agrippino), the shocks are identified up to a
#' Cholesky-type restriction within the block of `target` variables, so the order of `target`
#' matters. The impact matrix is normalised so that column `i` has impact `size[i]` on
#' `target[i]`. Missing instrument values must be zero-filled (this leaves the normalised
#' responses unchanged); NA are rejected.
#'
#' The extra output `shocks` (realised unit-variance structural shocks, \eqn{\varepsilon_t =
#' B_1'\Sigma_u^{-1}u_t}) and `reliability` (squared correlation of the instrument with the
#' shock; squared canonical correlations if `r > 1`) are sample objects computed from the
#' point estimates, valid when the instrument is relevant and exogenous; no censoring
#' adjustment is made.
#'
#' @param fit a [varm()] object.
#' @param m proxy: numeric vector, or `T x r` matrix of `r` proxies, aligned with
#'   the data (the first `p` rows are dropped to match the residuals).
#' @param target variable(s) (name or position) whose impact response is normalised to `size`;
#'   length `r`.
#' @param size size of the shock (impact response of `target`); scalar or length `r`.
#' @param horizon maximum horizon (impact is 0).
#' @param nboot number of bootstrap draws (`0` for point estimates only).
#' @param level confidence levels; bands are equal-tailed percentiles.
#' @param method `"mbb"` (default), `"wild"` or `"wild_normal"`; see Details.
#' @param block moving-block length for `method = "mbb"` (default Jentsch-Lunsford rule).
#' @return a `varm_proxy` object. For `r = 1`: `estimate` (`(H+1) x K`), `draws`
#'   (`(H+1) x K x nboot`), `bands` (`(H+1) x K x 2*length(level)`). For `r > 1` each gains a
#'   shock dimension after `K`. Also `gamma` (impact covariance with the proxy), `fstat` (see
#'   [proxy_fstat()]), `B1` (unit-variance impact matrix, `NULL` if not computable),
#'   `shocks`, `reliability`, `method`, `block`.
#' @references Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and corporate
#'   income tax changes in the United States. *AER* 103(4), 1212-1247.
#'
#'   Jentsch, C. and Lunsford, K. G. (2019). The dynamic effects of personal and corporate
#'   income tax changes in the United States: comment. *AER* 109(7), 2655-2678.
#'
#'   Jentsch, C. and Lunsford, K. G. (2022). Asymptotically valid bootstrap inference for
#'   proxy SVARs. *Journal of Business & Economic Statistics* (2022; working paper version: Cleveland Fed WP 16-19).
#' @seealso [varm_proxy_robust()], [varm_proxy_bayes()], [lp()] (LP-IV)
#' @examples
#' set.seed(1)
#' e <- matrix(rnorm(400), 200)
#' y <- matrix(0, 200, 2)
#' for (t in 2:200) y[t, ] <- 0.5 * y[t - 1, ] + c(e[t, 1] + 0.3 * e[t, 2], e[t, 2] + 0.5 * e[t, 1])
#' fit <- varm(y, p = 1)
#' varm_proxy(fit, m = e[, 1] + rnorm(200, sd = 0.5), horizon = 8, nboot = 100)
#' varm_proxy(fit, m = e[, 1] + rnorm(200, sd = 0.5), horizon = 8, nboot = 100, method = "mbb")
#' @export
varm_proxy <- function(
  fit,
  m,
  target = 1L,
  size = 1,
  horizon = 20L,
  nboot = 500L,
  level = c(0.68, 0.90),
  method = c("mbb", "wild", "wild_normal"),
  block = NULL
) {
  if (!inherits(fit, "varm")) {
    stop("`fit` must be a varm object", call. = FALSE)
  }
  method <- match.arg(method)
  p <- fit$p
  K <- fit$K
  H <- as.integer(horizon)
  M <- proxy_check_m(m, nrow(fit$Y))
  r <- ncol(M)
  iP <- resolve_var(target, fit$varnames, "target")
  if (length(iP) != r || anyDuplicated(iP)) {
    stop(
      "`target` must have one distinct variable per instrument",
      call. = FALSE
    )
  }
  size <- rep_len(size, r)
  Me <- M[-seq_len(p), , drop = FALSE]
  Tp <- nrow(fit$U)
  est <- proxy_irf(fit$B, fit$U, Me, iP, size, K, p, H)
  sn <- if (r == 1L) NULL else paste0("shock", seq_len(r))
  mom <- proxy_moments(fit$U, Me)
  B1 <- tryCatch(
    suppressWarnings(proxy_B1(mom$Su, mom$Sum, mom$Smm, iP)),
    error = function(e) NULL
  )
  if (!is.null(B1) && anyNA(B1)) {
    B1 <- NULL
  }
  out <- list(
    estimate = proxy_drop(est, r, list(0:H, fit$varnames, sn)),
    gamma = if (r == 1L) {
      setNames(mom$Sum[, 1], fit$varnames)
    } else {
      mom$Sum
    },
    horizon = H,
    p = p,
    K = K,
    r = r,
    varnames = fit$varnames,
    target = iP,
    size = size,
    level = level,
    fstat = proxy_fstat(fit, M, iP),
    B1 = B1,
    method = method
  )
  if (!is.null(B1)) {
    eps <- fit$U %*% solve(mom$Su, B1)
    out$shocks <- eps
    cm <- stats::cov(Me, eps)
    out$reliability <- if (r == 1L) {
      drop(stats::cor(Me, eps)^2)
    } else {
      sort(
        Re(
          eigen(solve(stats::cov(Me), cm %*% t(cm)), only.values = TRUE)$values
        ),
        decreasing = TRUE
      )
    }
  }
  if (nboot > 0L) {
    ell <- NA_integer_
    if (method == "mbb") {
      ell <- if (is.null(block)) round(5.03 * Tp^0.25) else as.integer(block)
      ell <- max(1L, min(ell, Tp - 1L))
    }
    out$block <- ell
    detp <- fit$det[-seq_len(p), , drop = FALSE]
    Y0 <- fit$Y[seq_len(p), , drop = FALSE]
    if (method == "mbb") {
      nb <- Tp - ell + 1L
      # recentring: mean of within-block position s over all nb blocks
      posmean <- function(A) {
        do.call(
          rbind,
          lapply(seq_len(ell), function(s) {
            colMeans(A[s:(s + nb - 1L), , drop = FALSE])
          })
        )
      }
      cU <- posmean(fit$U)
      cM <- posmean(Me)
      N <- ceiling(Tp / ell)
    }
    bd <- array(0, c(H + 1L, K, r, nboot))
    for (b in seq_len(nboot)) {
      if (method == "mbb") {
        st <- sample.int(nb, N, replace = TRUE)
        idx <- as.vector(outer(seq_len(ell), st, function(s, i) s + i - 1L))
        pos <- rep_len(seq_len(ell), N * ell)
        keep <- seq_len(Tp)
        Ub <- (fit$U[idx, , drop = FALSE] - cU[pos, , drop = FALSE])[
          keep,
          ,
          drop = FALSE
        ]
        Mb <- (Me[idx, , drop = FALSE] - cM[pos, , drop = FALSE])[
          keep,
          ,
          drop = FALSE
        ]
      } else {
        wb <- if (method == "wild") {
          sample(c(-1, 1), Tp, replace = TRUE)
        } else {
          stats::rnorm(Tp)
        }
        Ub <- fit$U * wb
        Mb <- Me * wb
      }
      Yb <- rbind(Y0, cpp_path(fit$B, K, p, Y0, detp, Ub))
      e <- cpp_var_ols(Yb, p, fit$det)
      bd[,,, b] <- proxy_irf(e$B, e$U, Mb, iP, size, K, p, H)
    }
    out <- proxy_bands(out, bd, level, fit$varnames, sn)
  }
  structure(out, class = "varm_proxy")
}

# ---- internals shared by the frequentist and Bayesian proxy routines --------

proxy_check_m <- function(m, Tn) {
  M <- as.matrix(m)
  storage.mode(M) <- "double"
  if (nrow(M) != Tn || anyNA(M)) {
    stop(
      "`m` must have no NAs and nrow(fit$Y) rows (zero-fill missing proxies)",
      call. = FALSE
    )
  }
  M
}

# uncentred second moments, as in Mertens-Ravn / Jentsch-Lunsford
proxy_moments <- function(U, M) {
  Tp <- nrow(U)
  list(
    Su = crossprod(U) / Tp,
    Sum = crossprod(U, M) / Tp,
    Smm = crossprod(M) / Tp
  )
}

# Mertens-Ravn (2013) B11/B21 identification, ported from the SVAR-IV MATLAB routine
# ProxySVARidentification.m: unit-variance impact matrix (K x r) of the r shocks.
proxy_B1 <- function(Su, Sum, Smm, iP) {
  K <- nrow(Su)
  r <- length(iP)
  perm <- c(iP, setdiff(seq_len(K), iP))
  S <- Su[perm, perm, drop = FALSE]
  Smu <- Sum[perm, , drop = FALSE]
  beta <- Smu %*% solve(Smm) # row i: regression of u_i on the instruments
  i1 <- seq_len(r)
  S11 <- S[i1, i1, drop = FALSE]
  S21 <- S[-i1, i1, drop = FALSE]
  S22 <- S[-i1, -i1, drop = FALSE]
  B21B11 <- beta[-i1, , drop = FALSE] %*% solve(beta[i1, , drop = FALSE])
  A <- S21 - B21B11 %*% S11
  Zeta <- B21B11 %*%
    S11 %*%
    t(B21B11) -
    (S21 %*% t(B21B11) + B21B11 %*% t(S21)) +
    S22
  B12B12 <- t(A) %*% solve(Zeta) %*% A
  B11B11 <- S11 - B12B12
  if (r == 1L) {
    B <- rbind(sqrt(B11B11[1, 1]), sqrt(B11B11[1, 1]) * B21B11)
  } else {
    B22B22 <- S22 + B21B11 %*% (B12B12 - S11) %*% t(B21B11)
    B12B22 <- (B12B12 %*% t(B21B11) + t(A)) %*% solve(B22B22)
    B11S1 <- diag(r) - B12B22 %*% B21B11
    B21S1 <- B21B11 %*% solve(B11S1)
    S1 <- t(chol(B11S1 %*% B11B11 %*% t(B11S1)))
    B <- rbind(solve(B11S1), B21S1) %*% S1
  }
  out <- matrix(0, K, r)
  out[perm, ] <- B
  out
}

# impact matrix (K x r), column i scaled to size[i] on variable iP[i]
proxy_impact <- function(Su, Sum, Smm, iP, size) {
  r <- length(iP)
  if (r == 1L) {
    # direction of Gamma only: no square root, so robust to noisy Sigma
    g <- Sum[, 1]
    return(matrix(size * g / g[iP], ncol = 1L))
  }
  B1 <- proxy_B1(Su, Sum, Smm, iP)
  sweep(B1, 2, size / B1[cbind(iP, seq_len(r))], "*")
}

# (H+1) x K x r responses for given reduced-form (B, moments)
proxy_irf_mom <- function(B, mom, iP, size, K, p, H) {
  G <- proxy_impact(mom$Su, mom$Sum, mom$Smm, iP, size)
  Th <- cpp_ma(B, K, p, H)
  out <- array(0, c(H + 1L, K, length(iP)))
  for (h in seq_len(H + 1L)) {
    out[h, , ] <- Th[,, h] %*% G
  }
  out
}

proxy_irf <- function(B, U, M, iP, size, K, p, H) {
  proxy_irf_mom(B, proxy_moments(U, M), iP, size, K, p, H)
}

proxy_drop <- function(a, r, dn) {
  if (r == 1L) {
    matrix(a, dim(a)[1], dim(a)[2], dimnames = dn[1:2])
  } else {
    dimnames(a) <- dn
    a
  }
}

# attach draws (H+1, K, r, N) and equal-tailed percentile bands to a varm_proxy list
proxy_bands <- function(out, bd, level, varnames, sn) {
  d <- dim(bd)
  r <- d[3]
  probs <- level_probs(level)
  bands <- aperm(apply(bd, 1:3, stats::quantile, probs = probs), c(2, 3, 4, 1))
  if (r == 1L) {
    out$draws <- array(bd, d[-3])
    out$bands <- array(bands, c(d[1:2], length(probs)))
    dimnames(out$bands) <- list(0:(d[1] - 1L), varnames, format(probs))
  } else {
    out$draws <- bd
    out$bands <- bands
    dimnames(out$bands) <- list(0:(d[1] - 1L), varnames, sn, format(probs))
  }
  out
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
#' @return named numeric `c(F = , F_robust = , F_hac = )`. With several instruments (`m` a matrix,
#'   `target` of the same length) a matrix with one row per target; `F`/`F_robust` are then joint
#' statistics over the instruments and `F_hac` is `NA` (the MOSW sets are single-instrument).
#' @seealso [varm_proxy_robust()]
#' @export
proxy_fstat <- function(fit, m, target = 1L, nwlags = NULL) {
  j <- resolve_var(target, fit$varnames, "target")
  M <- as.matrix(m)[-seq_len(fit$p), , drop = FALSE]
  one <- function(jj) {
    mod <- stats::lm(fit$U[, jj] ~ M)
    f <- unname(summary(mod)$fstatistic[1])
    fr <- if (requireNamespace("sandwich", quietly = TRUE)) {
      b <- stats::coef(mod)[-1]
      V <- sandwich::vcovHC(mod, type = "HC1")[-1, -1, drop = FALSE]
      unname(drop(crossprod(b, solve(V, b))) / length(b))
    } else {
      NA_real_
    }
    fh <- NA_real_
    if (ncol(M) == 1L) {
      cv <- proxy_cov(fit, M[, 1], nwlags)
      fh <- cv$Tp * cv$Gamma[jj]^2 / cv$W[cv$iG[jj], cv$iG[jj]]
    }
    c(F = f, F_robust = fr, F_hac = unname(fh))
  }
  res <- t(vapply(j, one, numeric(3)))
  if (length(j) == 1L) res[1, ] else `rownames<-`(res, fit$varnames[j])
}

#' @export
print.varm_proxy <- function(x, ...) {
  nm <- c(
    wild = "wild-bootstrap",
    wild_normal = "normal-wild-bootstrap",
    mbb = "moving-block-bootstrap",
    bayes = "posterior",
    arw = "posterior"
  )[[x$method]]
  cat(
    "<varm_proxy> shock to ",
    paste(x$varnames[x$target], collapse = ", "),
    " (size ",
    paste(x$size, collapse = ", "),
    "), horizon 0..",
    x$horizon,
    if (!is.null(x$draws)) {
      paste0(", ", dim(x$draws)[length(dim(x$draws))], " ", nm, " draws")
    },
    if (isTRUE(x$block > 0)) paste0(" (block length ", x$block, ")"),
    "\n",
    sep = ""
  )
  if (is.null(x$fstat)) {
    NULL
  } else if (is.matrix(x$fstat)) {
    cat("first-stage F by target:\n")
    print(signif(x$fstat[, 1:2, drop = FALSE], 4))
  } else {
    cat(
      "first-stage F = ",
      format(x$fstat[["F"]], digits = 4),
      "; robust F = ",
      format(x$fstat[["F_robust"]], digits = 4),
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$reliability)) {
    cat(
      "instrument reliability: ",
      paste(format(x$reliability, digits = 3), collapse = ", "),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}
