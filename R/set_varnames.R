#' @export
set_varnames <- function(object, nms, ...) {
  UseMethod("set_varnames")
}


#' @export
set_varnames.irf_varm <- function(object, nms, ...) {
  if(is.null(nms)) {
    return(object)
  }
  dnames_irf <- dimnames(object[[1]])
  dimnames(object$irfs) <- list(nms, nms, dnames_irf[[3]])

  dnames_birf <- dimnames(object[[2]])
  dimnames(object$boot_irfs) <- list(nms, nms, dnames_birf[[3]], dnames_birf[[4]])
  object
}
