
# mvtnorm::rmvnorm()
rmnorm <- function(n = 1, mean = rep(0, d), varcov, sqrt = NULL) {
  sqrt.varcov <- if (is.null(sqrt))
    chol(varcov)
  else sqrt
  d <- if (is.matrix(sqrt.varcov))
    ncol(sqrt.varcov)
  else 1
  mean <- outer(rep(1, n), as.vector(matrix(mean, d)))
  drop(mean + t(matrix(rnorm(n * d), d, n)) %*% sqrt.varcov)
}



#' Bootstrap `varm` model
#'
#'
boot_var <- function(obj, boot_fn = boot_resample) {
  resids <- base::scale(residuals(obj))
  boot_resids <- boot_fn(resids)
  sim_var(coefficients(obj), innov = boot_resids, nrep = nrow(boot_resids))
}



#' Simulate VAR
#'
#'
#' @references \url{https://github.com/MatthieuStigler/tsDyn/blob/a6aa2ebc041e448b3eb8598ffd117a83cfe316e5/tsDyn/R/VAR.sim.R}
#'
#' @seealso mlvar::simulateVAR  tsDyn::VAR.sim
#' @export
#' @examples
#'
#' obj <- varm(varm::spec(econdata::sw2001[,-1], .endo_lags = 4))
#' Acoef <- coefficients(obj)
#' set.seed(123)
#' dyn <- tsDyn:::VAR.sim(t(Acoef), n = 200, lag = 4, include = "const")
#' plot.ts(dyn)
#' set.seed(123)
#' ab <- gen_var(Acoef, nrep = 200)
#' plot.ts(ab)
#' waldo::compare(dyn, ab)
#'
#'
sim_var <- function(Acoef, innov = NULL, varcov = NULL, nrep = 200) {

  has_intercept <- attr(Acoef, "intercept")
  if(has_intercept) {
    B <- t(Acoef[-1,])
    ct <- Acoef[1,]
  }else{
    B <- t(Acoef)
    ct <- c(0,0,0)
  }
  K <- attr(Acoef, "K")
  p <- attr(Acoef, "p")

  varcov <- varcov %||% diag(1, nrow(B))
  innov <- innov %||% rmnorm(nrep, varcov = varcov)
  resb <- rbind(matrix(0, nrow = p, ncol = K), innov)

  simy <- matrix(0, ncol = K, nrow = nrep + p)
  for (i in (p + 1):(nrep + p)) {
    Y <- matrix(t(simy[i - c(1:p),, drop = FALSE]), ncol = 1)
    simy[i,] <- ct + B %*% Y + resb[i,]
  }
  simy[-c(1:p), , drop = FALSE]
}
