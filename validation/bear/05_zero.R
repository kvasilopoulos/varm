# BEAR irfres.m / qzerores.m zero restrictions vs varm id_zerosign.
# Deterministic check: a fully recursive set of impact zeros (shock2: 0 on var1; shock3: 0 on var1,var2) + positive diagonal must reproduce the Cholesky
# impact matrix exactly (both for BEAR's null-space algorithm replicated here and for varm) at the OLS point estimate.
source("validation/bear/common.R"); set.seed(5)
p <- 4; d <- load_bear_data(p); Y <- d$Y; K <- 3
f <- varm(Y, p); Dl <- t(chol(f$Sigma))
Z <- matrix(0, K, K); Z[1, 2] <- 1; Z[1:2, 3] <- 1          # [variable, shock]
S <- diag(K)
id <- id_zerosign(sign = S, zero = Z, nrot = 5, maxtry = 200)
v <- irf(f, id = id, horizon = 3, ci = "none")
a <- v$draws; cat("varm accepted rotations:", dim(a)[4], "\n")
mx <- max(sapply(seq_len(dim(a)[4]), function(r) max(abs(a[, , 1, r] - Dl))))
rep_("T12 varm zero+sign (recursive set) == Cholesky impact at OLS estimate", mx < 1e-9, "maxdiff", format(mx, digits = 3))
# BEAR qzerores replica (columns in decreasing # of zeros is ARRW; BEAR processes shocks 1..n in order and requires null-space) -> unique up to sign
bear_q <- function(Z, imp) { Qj <- NULL
  for (j in 1:K) { Zj <- if (any(Z[, j] == 1)) imp[Z[, j] == 1, , drop = FALSE] else NULL; Rj <- rbind(Zj, if (!is.null(Qj)) t(Qj))
    x <- rnorm(K); if (is.null(Rj)) q <- x / sqrt(sum(x^2)) else { N <- MASS::Null(t(Rj)); y <- N %*% (t(N) %*% x); q <- y / sqrt(sum(y^2)) }
    Qj <- cbind(Qj, q) }
  Qj }
# BEAR column order here: shock1 has no zeros, so it is drawn FIRST in the null space of nothing => random; shocks 2,3 then constrained.
# For a recursive scheme BEAR's sequential scheme needs the most-restricted shock first (a known order dependence of qzerores.m): test both orders.
Qb <- bear_q(Z, Dl)
cat("BEAR-replica (qzerores order 1,2,3, unrestricted shock first): rank-deficient/no admissible null space for shock 3:", anyNA(Qb) || length(Qb) == 0, "\n")
Zr <- Z[, 3:1]; Dq <- Dl %*% bear_q(Zr, Dl); Dq <- Dq[, 3:1]
sgn <- sign(diag(Dq)); Dq <- sweep(Dq, 2, sgn, `*`)
cat("BEAR-replica with most-restricted shock first: maxdiff to Cholesky:", format(max(abs(Dq - Dl)), digits = 3), "\n")
cat("NOTE: BEAR qzerores.m draws columns in shock order 1..n; if an unrestricted shock comes first the later columns can end up with no admissible null space (K-j < #zeros): ARRW (and varm) sort shocks by decreasing number of zeros.\n")
rep_("T12b varm accepts the same restrictions regardless of shock order (column permutation)", {
  perm <- c(3, 1, 2); idp <- id_zerosign(sign = S[perm, perm], zero = Z[perm, perm], nrot = 3, maxtry = 200)
  ap <- irf(f, id = idp, horizon = 0, ci = "none")$draws; ok <- all(sapply(seq_len(dim(ap)[4]), function(r) all(abs(ap[, , 1, r][Z[perm, perm] == 1]) < 1e-9))); ok })
