#' Least Square Estimation of VAR models
#'
#' @param x a `spec` object
#' @importFrom dplyr slice
#' @export
var_ls <- function(x) {

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
  resids <- Y - Z %*% A
  vcov <- crossprod(resids) / (N - K * p - p - 1)

  se <- crossprod(Z) %>%
    inv() %>%
    kronecker(vcov) %>%
    diag() %>%
    sqrt() %>%
    matrix(ncol = K) %>%
    `dimnames<-`(dimnames(coefs))

  tstat <- coefs/se
  pval <- 2*(1 - pt(tstat, df = N - K))

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
      init_data = x$lhs,
      K = K,
      p = p,
      N = N
    ),
    endo_varnames = x$endo_varnames,
    intercept = x$intercept,
    class = "var_ls"
  )
}


fast_ls <- function(x, p = 2) {
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



stable <- function(x) {
  Acomp <- comp(x$coefficients)
  roots <- eigen(Acomp)$values
  mod_roots <- Mod(roots)
  list(eigen = roots, mod_eigen = mod_roots)
}

is_stable <- function(x) {
  model <- stable(x)
  if (any(model$mod_roots < 1))
    FALSE
  else TRUE
}

check_stable <- function(x) {
  if (!is_stable(x)) {
    warning("model is not stationary")
  }
}


# methods -----------------------------------------------------------------



coef.var_ls <- function(object, ...) {
  structure(
    object$coefficients,
    intercept = attr(object, "intercept"),
    K = object$K,
    p = object$p
  )
}

print.var_ls <- function(x, digits = NULL) {
  cli::cat_rule(right = "var_ls()")
  cli::cat_line()
  coefx <- coef(x) %>%
    strip_attributes("K", "p", "intercept")
  format(coefx, digits = digits) %>%
    print(quote = FALSE)
}


