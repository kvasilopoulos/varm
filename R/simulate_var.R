
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

boot_var <- function(obj) {
  resids <- scale(obj$residuals)
  boot_resids <- boot_inst(resids)
  gen_var(obj$coefficients, innov = boot_resids, nrep = nrow(boot_resids))
}

#' Acoef <- var_ls(abvar::spec(econdata::sw2001[,-1], .endo_lags = 4))$coefficients
#' tsDyn::VAR.sim(t(Acoef[-1,]), n = 200, lag = 4, include = "none")
#' https://github.com/MatthieuStigler/tsDyn/blob/a6aa2ebc041e448b3eb8598ffd117a83cfe316e5/tsDyn/R/VAR.sim.R
gen_var <- function(Acoef, innov = NULL, varcov = NULL, nrep = 100) {
  B <- t(Acoef[-1,])
  K <- nrow(B)
  p <- ncol(B)/K

  varcov <- varcov %||% diag(1, nrow(B))
  innov <- innov %||% rmnorm(nsim, varcov = varcov)
  resb <- rbind(matrix(0, nrow = p, ncol = K), innov)

  simy <- matrix(0, ncol = K, nrow = nrep + p)
  for (i in (p + 1):(nrep + p)) {
    auxy <- matrix(t(simy[i - c(1:p),, drop = FALSE]), ncol = 1)
    simy[i,] <- B %*% auxy + resb[i,]
  }
  simy[-c(1:p), , drop = FALSE]
}


