#' Estimate a reduced-form VAR by least squares
#'
#' Equation-by-equation OLS (identical to Gaussian ML for the coefficients) of
#' \deqn{y_t = c + \sum_{l=1}^p A_l y_{t-l} + \Gamma w_t + u_t,\quad u_t \sim (0,\Sigma).}
#' The estimation loop is compiled (RcppArmadillo).
#'
#' @param Y data frame (preferred), `ts` or numeric matrix (T x K) of endogenous
#'   variables. Non-numeric columns of a data frame are dropped; a `Date` column
#'   is kept as `$dates` and used by the plots.
#' @param p lag order.
#' @param const include a constant.
#' @param trend include a linear trend.
#' @param exog optional T x q matrix or data frame of exogenous regressors
#'   (contemporaneous, treated like deterministic terms; future values are then
#'   required for forecasting).
#' @param vars <[`tidy-select`][tidyselect::language]> columns of `Y` that are
#'   the endogenous variables, e.g. `vars = c(infl, rate)`. Default: all numeric
#'   columns.
#' @return an object of class `varm` with elements `B` (m x K coefficient matrix,
#'   lags first then deterministics), `Sigma` (degrees-of-freedom corrected),
#'   `U` residuals, `XtXinv`, `logLik` and metadata. See [tidy.varm()],
#'   [glance.varm()] and [augment.varm()] for tibble summaries.
#' @details See `vignette("var-theory")` for definitions and the economic reading
#'   of the estimates.
#' @family estimation
#' @examples
#' fit <- varm(yraw, p = 2)
#' fit
#' summary(fit)
#' tidy(fit)
#' glance(fit)
#'
#' # choose the endogenous variables with tidyselect
#' varm(yraw, p = 1, vars = c(infl, rate))
#' @export
varm <- function(
  Y,
  p = 1L,
  const = TRUE,
  trend = FALSE,
  exog = NULL,
  vars = NULL
) {
  cl <- match.call()
  tsp_y <- if (stats::is.ts(Y)) stats::tsp(Y)
  Y <- as_data_matrix(Y, vars = rlang::enquo(vars))
  dates <- attr(Y, "dates")
  attr(Y, "dates") <- NULL
  p <- as.integer(p)
  if (p < 1L) {
    rlang::abort("`p` must be >= 1")
  }
  K <- ncol(Y)
  det <- make_det(nrow(Y), const, trend, exog)
  if (nrow(Y) - p <= K * p + ncol(det)) {
    rlang::abort(paste0("not enough observations for p = ", p))
  }
  est <- cpp_var_ols(Y, p, det)
  vn <- colnames(Y)
  rn <- coef_names(vn, p, colnames(det))
  dimnames(est$B) <- list(rn, vn)
  dimnames(est$Sigma) <- dimnames(est$Sigma_ml) <- list(vn, vn)
  colnames(est$U) <- vn
  dimnames(est$XtXinv) <- list(rn, rn)
  fields <- est[c("B", "Sigma", "Sigma_ml", "U", "XtXinv", "logLik")]
  structure(
    c(
      fields,
      list(
        Y = Y,
        det = det,
        p = p,
        K = K,
        nobs = est$nobs,
        npar = est$npar,
        varnames = vn,
        detnames = colnames(det),
        dates = dates,
        tsp = tsp_y,
        call = cl
      )
    ),
    class = "varm"
  )
}

#' Select the VAR lag order
#'
#' Computes AIC, HQ (Hannan-Quinn), SC (Schwarz/BIC) and FPE for lags `1..pmax`
#' on a common estimation sample (the first `pmax` observations are dropped for
#' every `p`).
#' @inheritParams varm
#' @param pmax maximum lag.
#' @return list with `selection` (the optimal lag per criterion) and the
#'   `criteria` matrix.
#' @family estimation
#' @examples
#' var_select(yraw, pmax = 4)$selection
#' @export
var_select <- function(
  Y,
  pmax = 8L,
  const = TRUE,
  trend = FALSE,
  exog = NULL,
  vars = NULL
) {
  Y <- as_data_matrix(Y, vars = rlang::enquo(vars))
  attr(Y, "dates") <- NULL
  det <- make_det(nrow(Y), const, trend, exog)
  ic <- cpp_var_ic(Y, as.integer(pmax), det)
  dimnames(ic) <- list(seq_len(pmax), c("AIC", "HQ", "SC", "FPE"))
  list(selection = apply(ic, 2, which.min), criteria = ic)
}

# accessors -----------------------------------------------------------------

#' @export
coef.varm <- function(object, ...) object$B

#' @export
resid.varm <- function(object, ...) object$U
#' @export
residuals.varm <- resid.varm

#' @export
fitted.varm <- function(object, ...) {
  object$Y[-seq_len(object$p), , drop = FALSE] - object$U
}

#' @export
vcov.varm <- function(object, ...) object$Sigma

#' @export
nobs.varm <- function(object, ...) object$nobs

#' @export
logLik.varm <- function(object, ...) {
  structure(
    object$logLik,
    df = object$npar * object$K + object$K * (object$K + 1) / 2,
    nobs = object$nobs,
    class = "logLik"
  )
}

#' Lag coefficient matrices
#'
#' @param object a `varm`/`bvarm` object or a coefficient matrix.
#' @param ... unused.
#' @return array `[K, K, p]` with `A[,,l]` the coefficient matrix on lag `l`
#'   (row = equation, column = regressor variable).
#' @family estimation
#' @examples
#' lag_matrices(varm(yraw, p = 2))[, , 1]
#' @export
lag_matrices <- function(object, ...) UseMethod("lag_matrices")
#' @export
lag_matrices.varm <- function(object, ...) {
  lag_matrices_B(object$B, object$K, object$p, object$varnames)
}
lag_matrices_B <- function(B, K, p, vn = NULL) {
  A <- array(0, c(K, K, p), dimnames = list(vn, vn, NULL))
  for (l in seq_len(p)) {
    A[,, l] <- t(B[(l - 1) * K + seq_len(K), , drop = FALSE])
  }
  A
}

#' @export
print.varm <- function(x, digits = 4, ...) {
  cat(
    "<varm> VAR(",
    x$p,
    ") with ",
    x$K,
    " variables, ",
    x$nobs,
    " effective observations\n",
    sep = ""
  )
  cat(
    "logLik = ",
    format(x$logLik, digits = 6),
    "; max |eigenvalue| = ",
    format(max(stability(x)), digits = 4),
    "\n\n",
    sep = ""
  )
  print(round(x$B, digits))
  invisible(x)
}

#' @export
summary.varm <- function(object, ...) {
  se <- sqrt(outer(diag(object$XtXinv), diag(object$Sigma)))
  tval <- object$B / se
  df <- object$nobs - object$npar
  pval <- 2 * stats::pt(-abs(tval), df)
  Yd <- object$Y[-seq_len(object$p), , drop = FALSE]
  rss <- colSums(object$U^2)
  tss <- colSums(sweep(Yd, 2, colMeans(Yd))^2)
  structure(
    list(
      object = object,
      se = se,
      tvalue = tval,
      pvalue = pval,
      r.squared = 1 - rss / tss,
      adj.r.squared = 1 - (rss / df) / (tss / (object$nobs - 1))
    ),
    class = "summary.varm"
  )
}

#' @export
print.summary.varm <- function(x, digits = 4, ...) {
  o <- x$object
  cat("<varm> VAR(", o$p, "), K = ", o$K, ", T = ", o$nobs, "\n", sep = "")
  for (k in seq_len(o$K)) {
    cat(
      "\nEquation: ",
      o$varnames[k],
      "   R2 = ",
      format(x$r.squared[k], digits = 4),
      "  adj.R2 = ",
      format(x$adj.r.squared[k], digits = 4),
      "\n",
      sep = ""
    )
    tab <- cbind(
      Estimate = o$B[, k],
      `Std. Error` = x$se[, k],
      `t value` = x$tvalue[, k],
      `Pr(>|t|)` = x$pvalue[, k]
    )
    stats::printCoefmat(tab, digits = digits, has.Pvalue = TRUE)
  }
  cat("\nResidual covariance:\n")
  print(round(o$Sigma, digits + 2))
  cat("Residual correlation:\n")
  print(round(stats::cov2cor(o$Sigma), 3))
  invisible(x)
}

#' Eigenvalue moduli of the companion matrix
#'
#' A VAR is covariance stationary iff all moduli are < 1.
#' @param object a `varm` or `bvarm` object.
#' @return numeric vector of moduli (descending); for Bayesian objects, a matrix `draws x Kp`.
#' @family estimation
#' @examples
#' stability(varm(yraw, p = 2))
#' @export
stability <- function(object) UseMethod("stability")
#' @export
stability.varm <- function(object) cpp_stability(object$B, object$K, object$p)
