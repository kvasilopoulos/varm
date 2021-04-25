#' Least Square Estimation of VAR models
#'
#' @param x a `spec` object
#' @importFrom dplyr slice
#' @export
varm <- function(x) {

  # TODO no need to store so much metadata

  p <- x$endo_lags
  K <- NCOL(x$lhs)
  N <- NROW(x$lhs)

  has_trend <- x$spec_intercept
  has_intercept <- x$spec_intercept

  Y <- x$lhs[-c(1:p), ] %>%
    as.matrix()
  Z <- x$rhs[-c(1:p), ] %>%
    as.matrix()

  coefs <- A <- inv(crossprod(Z)) %*% crossprod(Z, Y)
  fitted <- Z %*% A

  resids <- Y - fitted
  vcov <- crossprod(resids) / (N - K * p - p - 1)

  # for summary
  rss <- colSums(resids^2)
  mss <- colSums(fitted^2)
  rsquared <- mss / (rss + mss)
  # TODO df <-
  # TODO Fstat
  # TODO LogLikk
  # TODO se equation
  # TODO adj_rsquared <- 1 - (1 - rsquared) * n / df

  se <- crossprod(Z) %>%
    inv() %>%
    kronecker(vcov) %>%
    diag() %>%
    sqrt() %>%
    matrix(ncol = K) %>%
    `dimnames<-`(dimnames(coefs))

  tstat <- coefs/se
  pval <- 2*pt(-abs(tstat), df = N - K)

  # confidence intervals
  tcrit <- qt(0.95, N)
  conf  <- list(lower_ci = coefs - tcrit*se, upper_ci =  coefs + tcrit*se)

  structure(
    list(
      coefficients = meta_matrix(coefs, intercept = x$intercept, K = K, p = p),
      residuals = resids,
      lhs = x$lhs,
      vcov = vcov,
      se = se,
      tstat = tstat,
      pvalue = pval,
      confint = conf,
      r_squared = rsquared,
      init_data = x$lhs,
      K = K,
      p = p,
      N = N
    ),
    endo_varnames = x$endo_varnames,
    intercept = x$intercept,
    class = "varm"
  )
}

#' Fast least square estimation
#'
#' Stores fewer metadata
#'
#' @export
varm_fast <- function(x, p = 2) {
  matx <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  Z <- mat_lag(matx, lags = p)
  Y <- matx[-c(1:p),]
  A <- inv(crossprod(Z)) %*% crossprod(Z, Y)
  U <- Y - Z %*% A
  vcov <- crossprod(U) / (nr - nc * p - p - 1)
  list(
    coefficients = A,
    residuals = U,
    vcov = vcov,
    K = K,
    p = p
  )
}




# methods -----------------------------------------------------------------

coef.varm <- function(object, ...) {
  structure(
    object$coefficients,
    intercept = attr(object, "intercept"),
    K = object$K,
    p = object$p
  )
}

#' @export
print.varm <- function(x, digits = getOption("digits") - 3, ...) {
  cli::cat_rule(left = "varm()")
  cli::cat_line()
  coefx <- coef(x) %>%
    strip_attributes()
  format(coefx, digits = digits) %>%
    print(quote = FALSE, ...)
}


#' @export
#' @importFrom stats coef coefficients printCoefmat pt qt quantile resid residuals runif sd
summary.varm <- function(x, digits = max(3L, getOption("digits") - 3L),
                           signif.stars = getOption("show.signif.stars"),
                           signif.legend = signif.stars, ...) {
  coefx <- coef(x) %>%
    strip_attributes()
  nms <- colnames(coefx)
  se <- x$se
  tstat <- x$tstat
  pval <- x$pvalue
  ncoef <- ncol(coefx)
  signif_legend <- c(rep(FALSE, ncoef - 1), signif.legend = signif.legend)
  for (i in 1:ncoef) {
    cli::cat_rule(left = nms[i])
    coefs <- cbind(
      Estimate = coefx[, i],
      `Std. Error` = se[, i],
      `t value` = tstat[, i],
      `Pr(>|t|)` = pval[, i]
    )
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               signif.legend = signif_legend[i])
  }
  cli::cat_line()
  cli::cat_rule(left = "sigma")
  print(x$vcov, digits = digits)
}
