# Internal helpers -----------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

DATE_CLASSES <- c("Date", "POSIXt", "yearmon", "yearqtr")

is_date_col <- function(x) inherits(x, DATE_CLASSES)

# Coerce data (data frame, ts, matrix, vector) to a numeric matrix. Non-numeric
# columns of a data frame are dropped; a date column is kept as attribute
# "dates". `vars` is an optional quosure (tidyselect) choosing the columns.
as_data_matrix <- function(Y, arg = "Y", vars = NULL) {
  label <- deparse(substitute(Y))[1]
  dates <- NULL
  if (is.data.frame(Y)) {
    date_col <- Filter(is_date_col, Y)
    dates <- if (length(date_col)) date_col[[1]]
    if (inherits(dates, c("yearmon", "yearqtr"))) {
      dates <- as.Date(dates)
    }
    Y <- df_to_matrix(Y, vars, arg)
  } else {
    if (is.null(dim(Y))) {
      Y <- matrix(Y, ncol = 1L, dimnames = list(NULL, label))
    }
    Y <- select_columns(Y, vars, arg)
  }
  if (!is.numeric(Y)) {
    rlang::abort(paste0("`", arg, "` must be numeric"))
  }
  if (anyNA(Y)) {
    rlang::abort(paste0("`", arg, "` contains missing values"))
  }
  storage.mode(Y) <- "double"
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  }
  attr(Y, "dates") <- dates
  Y
}

select_columns <- function(Y, vars, arg) {
  if (is.null(vars) || rlang::quo_is_null(vars)) {
    return(Y)
  }
  cn <- colnames(Y) %||% paste0("y", seq_len(ncol(Y)))
  Y[, select_vars(vars, cn, "variable"), drop = FALSE]
}

df_to_matrix <- function(df, vars, arg) {
  if (is.null(vars) || rlang::quo_is_null(vars)) {
    keep <- vapply(df, function(x) is.numeric(x) && !is_date_col(x), NA)
    df <- df[keep]
  } else {
    df <- df[select_vars(vars, names(df), "variable")]
  }
  if (ncol(df) == 0L) {
    rlang::abort(paste0("`", arg, "` has no numeric columns"))
  }
  as.matrix(df)
}

# Deterministic regressors (T x d): constant, trend, exogenous columns
make_det <- function(n, const, trend, exog) {
  d <- matrix(0, n, 0L)
  if (const) {
    d <- cbind(d, const = rep(1, n))
  }
  if (trend) {
    d <- cbind(d, trend = seq_len(n))
  }
  if (!is.null(exog)) {
    exog <- as_data_matrix(exog, "exog")
    if (nrow(exog) != n) {
      rlang::abort("`exog` must have the same number of rows as `Y`")
    }
    if (!all(nzchar(colnames(exog)))) {
      colnames(exog) <- paste0("exog", seq_len(ncol(exog)))
    }
    d <- cbind(d, exog)
  }
  d
}

coef_names <- function(varnames, p, detnames) {
  c(unlist(lapply(seq_len(p), function(l) paste0(varnames, ".l", l))), detnames)
}

# Backward-compatible resolver (positions from names/positions/tidyselect).
resolve_var <- function(x, varnames, arg = "variable") {
  select_vars(x, varnames, arg)
}

# Quantile probabilities for central intervals of the given levels
level_probs <- function(level) {
  sort(unique(c((1 - level) / 2, 1 - (1 - level) / 2)))
}

# Pair up the sorted probabilities: widest interval first.
# Returns the level of each (lower, upper) pair of `probs`.
band_levels <- function(probs) {
  n <- length(probs) / 2
  round(probs[length(probs) + 1 - seq_len(n)] - probs[seq_len(n)], 10)
}

# Quantiles of forecast draws f [N, h, K] over the draws: array [h, K, probs]
forecast_bands <- function(f, probs, varnames) {
  b <- aperm(apply(f, c(2, 3), stats::quantile, probs = probs), c(2, 3, 1))
  dimnames(b) <- list(NULL, varnames, format(probs))
  b
}

# cpp output K x K x (H1 * N) -> 4-d array
as_irf4 <- function(cube, K, H1) {
  array(cube, dim = c(K, K, H1, length(cube) / (K * K * H1)))
}
