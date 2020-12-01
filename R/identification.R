# TODO make them all methods

id_shock <- function(type = c("none", "cholesky", "triangular",
                              "sign", "magnitude", "zero")) {
  type <- match.arg(type)
  id_fun <-
    switch(
      type,
      none = id_none,
      choleski = id_chol,
      triangular = id_triangular,
      sign = id_sign,
      magnitude = id_magnitude,
      zero = id_zero)
}

id_none <- function(B = NULL,  ...) {
  function(B) {
    B
  }
}

id_chol <- function(B = NULL, ...) {
  function(B) t(chol(B, ...))
}


#' Identification with triangluar decomposition
#'
#' @examples
#' \dontrun{
#' A <- matrix(c(1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1), 4)
#' luDec <- lu(A)
#' L <- expand(luDec)$L
#' U <- expand(luDec)$U
#' p <- expand(luDec)$P
#' # LU is a row-permuted version of A
#' L %*% U
#' #Going back to the original identity A = PLU we can recover A
#' P %*% L %*% U
#' }
#'
#' @importFrom Matrix expand lu
#' @export
id_triangular <- function(x, ...) {
  function(x) {
    out <- t(as.matrix(Matrix::expand(Matrix::lu(x, ...))$U))
    dimnames(out) <- dimnames(x)
    out
  }
}

id_lr <- function() {
  # Finf_big = inv(eye(length(Fcomp))-Fcomp); % from the companion
  # Finf = Finf_big(1:nvar,1:nvar);
  # D  = chol(Finf*sigma*Finf')'; % identification: u2 has no effect on y1 in the long run
  #           invA = Finf\D;
}


# There is work to be done on rotation matrix
id_sign <- function(x, sign_vec) {
  function(B) {
  # out <- t(chol(B))
  # t(out) %*% t(sign_vec)
  }
}

id_magnitude <- function(x) {

}

id_zero <- function(x) {

}

id_iv <- function(x, iv) {

  # component <- match.arg(component)
  mres <- model$residuals
  t <- nrow(mres)
  n <- ncol(mres)
  p <- model$p
  bet <- model$coefficients

  matM <- as.matrix(iv)[-c(1:p),, drop = FALSE]

  # Identification
  gamma <- crossprod(mres, matM) / t
}
id_proxy <- id_iv
