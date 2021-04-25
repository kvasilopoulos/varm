bvar_minn <- function(Yraw, p = 4, nsave = 10000, nburn = 0) {

  nms <- colnames(Yraw)
  mats <- spec(Yraw, .endo_lags = p)
  ols_dat <- varm(mats)
  ntot <- nsave + nburn
  Traw <- NROW(Yraw)
  Ylag <- as.matrix(mats$lhs_lags)
  M <- NCOL(Yraw)

  # # Matrices ----------------------------------------------------------------

  Y <- as.matrix(mats$lhs)[-c(1:p), ]
  X <- as.matrix(mats$rhs)[-c(1:p), ]
  Z <- kronecker(eye(M), Y) # Block diagonal matrix Z

  K <- NCOL(X)
  TT <- Traw - p

  # # Priors ------------------------------------------------------------------

  A_OLS <- ols_dat$coefficients
  a_OLS <- c(A_OLS)
  SIGMA_OLS <- ols_dat$vcov

  # # Hyperparameters
  A_prior <- zeros(K, M)
  # diag(A_prior) <- 1 #if(rw == TRUE)
  a_prior <- c(A_prior)

  # # prior_hyper -------------------------------------------------------------

  a_bar_1 = 1 # own lags
  a_bar_2 = 1 # off diagnonal lags
  a_bar_3 = 10^6 # on exogenous variables

  sigma_sq <-  rep(0, M)
  alpha_i <- zeros(p, M)

  sigma_sq  <-  zeros(M,1) # % vector to store residual variances
  # Ylags <- ones(N - p,1) #;

  # Rearrange the X mat
  ind <- c(1, matrix(2:13, 4, 3, byrow = TRUE) %>% c())
  Ylags <- X[,ind]

  # alpha_i <- zeros(K, M)
  sigma_sq <- rep(0, M)
  for (i in 1:M) {
    # Dependent variable in i-th equation
    Y_i = Y[, i]
    # OLS estimates of i-th equation
    alpha_i = inv(t(Ylags) %*% Ylags) %*% (t(Ylags) %*% Y_i)
    sigma_sq[i] = (1/(Traw - p + 1)) %*% crossprod(Y_i - Ylags %*% alpha_i)
  }

  V_i <- zeros(M, M)
  V_pr <- a_bar_3 * ones(M, 1) #rep(0.9, 3)

  for (l in 1:p) {
    for (i in 1:M) {
      for (j in 1:M) {
        if (i == j) {
          V_i[i,j] <-  a_bar_1^2 / (l^2)
        }else{
          V_i[i,j] <-  (a_bar_1^2)*(a_bar_2^2)*sigma_sq[i] / ((l^2)*sigma_sq[j])
        }
      }
    }
    V_pr <- cbind(V_pr, V_i)
  }

  V_pr <- t(V_pr)
  V_prior <- diag(vec(V_pr))
  SIGMA <- diag(sigma_sq)

  # # Posteriors --------------------------------------------------------------

  # % Storage space for posterior draws
  alpha_draws <- zeros(nsave, K*M)
  ALPHA_draws <- array(0, dim = c(nsave, K, M))
  SIGMA_draws <- array(0 , dim = c(nsave, M, M))

  # if (irep > nburn) {
  #   Bv <- array(0, dim = c(M, M, p))
  #   for (i_1 in 1:p) {
  #     alpha_index <- (1 + (i_1 - 1) * M):(i_1 * M)
  #     alpha_index <- alpha_index + 1 # add one for constant
  #     Bv[, , i_1] <- ALPHA[alpha_index, ]
  #   }
  #   shock <- t(chol(SIGMA))
  #   d <- diag(diag(shock))
  #   shock <- solve(d) %*% shock
  #   all_responses[irep - nburn, , , ] <- impulse(Bv, shock, ihor)
  # }


  for (irep in 1:ntot ) {
    alpha <- NA
    for (i in 1:M) { # something wrong is here
      ksi <- inv(V_prior[((i - 1) * K + 1):(i * K), ((i - 1) * K + 1):(i * K)])
      V_post <- inv(ksi + t(X) %*% X / SIGMA[i, i])
      a_post <- V_post %*% (ksi %*% a_prior[((i - 1) * K + 1):(i * K)] + t(X) %*% Y[, i]/SIGMA[i, i])
      alpha[((i - 1) * K + 1):(i * K)] <- a_post + t(chol(V_post)) %*% stats::rnorm(K)
    }
    ALPHA <- matrix(alpha, K, M)
    if (irep > nburn) {
      alpha_draws[irep - nburn, ] <- alpha
      ALPHA_draws[irep - nburn, , ] <- ALPHA
      SIGMA_draws[irep - nburn, , ] <- SIGMA
    }
  }

  # Posterior mean of parameters:
  ALPHA_mean = apply(ALPHA_draws, c(2,3), mean) # posterior mean of ALPHA
  SIGMA_mean = apply(SIGMA_draws, c(2,3), mean) #posterior mean of SIGMA

  # %Posterior standard deviations of parameters:
  ALPHA_std = apply(ALPHA_draws, c(2,3), sd) # posterior mean of ALPHA
  SIGMA_std = apply(SIGMA_draws, c(2,3), sd) #posterior mean of SIGMA


  structure(
    list(
      Yraw = Yraw,
      SIGMA = SIGMA,
      ALPHA = ALPHA,
      ALPHA_draws = ALPHA_draws,
      SIGMA_draws = SIGMA_draws,
      M = M,
      p = p,
      nms = nms
    ),
    class = "bvar_minn"
  )
}


#' @export
irf.bvar_minn <- function(x, nstep = 20, ...) {
  # irf ---------------------------------------------------------------------

  M <- x$M
  p <- x$p
  ALPHA <- x$ALPHA
  SIGMA <- x$SIGMA
  nms <- x$nms

  # think of storing to a 4-dim BV and then squeeze outside the loop
  Bv = array(0, dim = c(M, M, p))
  # This is the partition function fix later ~ insted of list do 3-dim array
  for (i in 1:p) {
    Bv[,,i] = ALPHA[(1 + ((i - 1)*M + 1)):(i * M + 1), ];
  }
  shock <- chol(SIGMA) #

  irf_bvar <- function(By, smat, nstep) {
    dims <- dim(By)
    nr <- dims[1]
    nc <- dims[2]
    nlag <- dims[3]
    response <- array(0, c(nc, nr, nstep))
    response[ , , 1] <- t(smat) # no reason to transpose since it is diagonal
    for (it in 2:nstep) {
      for (ilag in 1:min(nlag, it - 1)) {
        response[, , it] = response[, , it] + By[, , ilag] %*% response[, , it - ilag]
      }
    }
    response
  }
  x <- irf_bvar(Bv, shock, 20)
  dimnames(x) <- list(nms, nms, 1:20)
  class(x) <- "irf_bvar_minn"
  x
}

#' @export
autoplot.irf_bvar_minn <- function(x, ...) {
  x %>%
    array_to_tbl() %>%
    ggplot(aes(horizon, irf)) +
    geom_line() +
    facet_wrap(response ~impulse, scales = "free_y") +
    theme_default()

}

