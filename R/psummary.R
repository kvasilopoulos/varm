#' @export
psummary <- function(x, ...) {
  UseMethod("psummary")
}


#' Pretty Summary
#'
#' @export
psummary.varm <- function(object, include = c("none", "tstat", "pvalue", "both"),
                     digits = max(3L, getOption("digits") - 3L),
                     ...) {

  z <- object
  z$cnames <- attr(object, "endo_varnames")
  nc <- z$K
  p <- z$p
  intercept <- attr(object, "intercept")
  num_int <- mean(intercept)

  include <- match.arg(include)
  # nwidth <- max(sapply(object$cnames, nchar))

  # we make 2 loop for preallocation
  if (include == "none") {
    nr <- z$K*z$p*2 + num_int*2
    slice <- seq(1, nr, 2)
  }else if (include %in% c("tstat","pvalue")) {
    nr <- z$K*z$p*3 + num_int*2
    slice <- seq(1, nr, 3)
  }else if (include == "both") {
    nr <- z$K*z$p*4 + num_int*2
    slice <- seq(1, nr, 4)
  }

  lagnames <- paste0(z$cnames, "(-", rep(1:p, each = nc) , ")")
  new_mat <- matrix(NA, nr, nc, dimnames = list(rep("", nr ), z$cnames))
  rownames(new_mat)[slice] <- c("Intercept", lagnames)

  if (include == "none") {
    new_mat[slice,] <- format(z$coefficients, digits = digits, width = 8, justify = "centre")
    new_mat[-slice,] <- paste0("(", format(z$se, digits = digits, justify = "centre"),")")
  }else if (include %in% c("tstat","pvalue")) {
    new_mat[seq(1, nr, 3),] <- format(z$coefficients, digits = digits,width = 8, justify = "centre")
    new_mat[seq(2, nr, 3),] <- paste0("(", format(z$se, digits = digits, width = 7, justify = "centre"),")")
    new_mat[seq(3, nr, 3),] <- paste0("[", format(
      if (include == "tstat") z$tstat else z$pvalue, digits = digits,
      width = 6, justify = "left"),"]")

  }else if (include == "both") {
    # num1 <- seq(1, z$K*z$p*2, 4)
    # num2 <- seq(2, z$K*z$p*2, 4)
    # num3 <- seq(3, z$K*z$p*2, 4)
    # num3 <- seq(4, z$K*z$p*2, 4)
    # nr <- z$K*z$p*4
  }
  # ret_mat <-
  ans <- list()
  ans$mat <- new_mat
  ans$nr <- nr
  ans$slice <- slice
  class(ans) <-  "psummary.varm"
  return(ans)
}

#' @export
print.psummary.varm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {

  cat("\nStandard errors in () & t-statistic in []\n\n")
  print.default(x$mat, digits = digits, print.gap = 2L, quote = FALSE, right = T)

}

