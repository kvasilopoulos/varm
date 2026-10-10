# Internal helpers ---------------------------------------------------------

as_data_matrix <- function(Y, arg = "Y") {
  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
  }
  if (is.null(dim(Y))) {
    Y <- matrix(Y, ncol = 1L, dimnames = list(NULL, deparse(substitute(Y))))
  }
  if (!is.numeric(Y)) {
    stop("`", arg, "` must be numeric", call. = FALSE)
  }
  if (anyNA(Y)) {
    stop("`", arg, "` contains missing values", call. = FALSE)
  }
  storage.mode(Y) <- "double"
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  }
  Y
}

# Deterministic regressors (T x d): constant, trend, exogenous columns
make_det <- function(n, const, trend, exog) {
  d <- NULL
  if (const) {
    d <- cbind(d, const = rep(1, n))
  }
  if (trend) {
    d <- cbind(d, trend = seq_len(n))
  }
  if (!is.null(exog)) {
    exog <- as_data_matrix(exog, "exog")
    if (nrow(exog) != n) {
      stop("`exog` must have the same number of rows as `Y`", call. = FALSE)
    }
    if (is.null(colnames(exog)) || !all(nzchar(colnames(exog)))) {
      colnames(exog) <- paste0("exog", seq_len(ncol(exog)))
    }
    d <- cbind(d, exog)
  }
  if (is.null(d)) {
    d <- matrix(0, n, 0L)
  }
  d
}

coef_names <- function(varnames, p, detnames) {
  c(unlist(lapply(seq_len(p), function(l) paste0(varnames, ".l", l))), detnames)
}

resolve_var <- function(x, varnames, arg = "variable") {
  if (is.character(x)) {
    i <- match(x, varnames)
    if (anyNA(i)) {
      stop(
        "unknown ",
        arg,
        ": ",
        paste(x[is.na(i)], collapse = ", "),
        call. = FALSE
      )
    }
    i
  } else {
    as.integer(x)
  }
}

# set 4-d dim on cpp output K x K x (H1*N)
as_irf4 <- function(cube, K, H1) {
  array(cube, dim = c(K, K, H1, length(cube) / (K * K * H1)))
}
