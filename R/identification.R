id_shock <- function(type = c("none", "choleski", "triangular",
                              "sign", "magnitude", "zero")) {
  type <- match.arg(type)
  id_fun <-
    switch(
      type,
      none = id_none(),
      choleski = id_chol,
      triangular = id_triangular,
      sign = id_sign,
      magnitude = id_magnitude,
      zero = id_zero)
}

id_none <- function(x = NULL, ...) {
  function(x) diag(nrow(x), ...)
}

id_chol <- function(x = NULL, ...) {
  function(x) t(chol(x, ...))
}


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
id_triangular <- function(x, ...) {
  function(x) Matrix::expand(Matrix::lu(x, ...))$L
}

id_sign <- function(x) {

}

id_magnitude <- function(x) {

}

id_zero <- function(x) {

}
