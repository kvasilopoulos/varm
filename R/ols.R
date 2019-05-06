#' Least Square Estimation of VAR models
#'
#' @param x a `spec` object
#' @importFrom dplyr slice
#' @export
ols <- function(x) {

  p <- x$endo_lags
  K <- NCOL(x$lhs)
  N <- NROW(x$lhs)

  has_trend <- x$spec_intercept
  has_intercept <- x$spec_intercept


  Y <- x$lhs %>% slice(-(1:p)) %>% as.matrix()
  Z <- x$rhs %>% slice(-(1:p)) %>% as.matrix()

  coefs <- A <- inv(crossprod(Z)) %*% crossprod(Z, Y)
  resids <- Y - Z %*% A
  vcov <- crossprod(resids) / (N - K * p - p - 1)

  structure(
    list(
      A = coefs,
      resids = resids,
      vcov = vcov,
      K = K,
      p = p,
      N = N
    ),
    class = "ols"
  )
}
