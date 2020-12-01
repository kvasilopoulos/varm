

#' This happens in id_iv and id_sign
#'
irf_partial <- function(A, B, K, p, h) {

  irs <- matrix(0, irhor + p, n, dimnames = list(NULL, colnames(model$matY)))
  irs[p + 1, ] <- gamma / gamma[1]
  for (jj in 2:irhor) {
    lvars <- c(t(irs[(p + jj - 1):jj, ])) # collapse to vector
    irs[p + jj, ] <- lvars %*% bet[-1, ] # remove constant
  }
  irs <- irs[-c(1:p), ]
  irs
}

# se <- function() {}
# se_boot <- function(engine, nboot = 100, seed = NULL) {}



#' Impulse Response Function Algorithm 2 (this may be wrong)
#'
irf_algo3 <- function(Acomp, B, K, p, h) {
  nms <- rownames(B)
  irf_comp <- array(0, c(K * p, h + 1, K * p))
  irf_comp[1:K, 1, 1:K] <- B
  for (i in 1:h) {
    irf_comp[, i + 1, ] <- irf_comp[, i, ] %*% t(Acomp)
  }
  irf_comp <- irf_comp[1:K, , 1:K]
  out <- aperm(irf_comp, c(1, 3, 2))
  dimnames(out) <- list(nms, nms, NULL)
  # return(list(out,B))
  out
}


create_nms <- function(var_nm) {
  len <- length(var_nm)
  nms <- matrix(NA, len, len)
  for (i in 1:len) {
    for (j in 1:len)
      nms[i, j] <- paste(var_nm[j], var_nm[i], sep = " ~ ")
  }
  c(nms)
}

irf_algo2 <- function(A, B, K, p, h) {

  if (h >= p) {
    As <- array(0, dim = c(K, K, nstep + 1))
    for (i in (p + 1):(nstep + 1)) {
      As[, , i] <- matrix(0, nrow = K, ncol = K)
    }
  }
  else {
    As <- array(0, dim = c(K, K, p))
  }
  for (i in 1:p) {
    As[, , i] <- A[[i]]
  }
  Phi <- array(0, dim = c(K, K, nstep + 1))
  Phi[, , 1] <- base::diag(K)
  Phi[, , 2] <- Phi[, , 1] %*% As[, , 1]
  if (nstep > 1) {
    for (i in 3:(nstep + 1)) {
      tmp1 <- Phi[, , 1] %*% As[, , i - 1]
      tmp2 <- matrix(0, nrow = K, ncol = K)
      idx <- (i - 2):1
      for (j in 1:(i - 2)) {
        tmp2 <- tmp2 + Phi[, , j + 1] %*% As[, , idx[j]]
      }
      Phi[, , i] <- tmp1 + tmp2
    }
  }
}



