partition <- function(x, transpose = TRUE, prefix_names = "lag",
                      keep_intercept = TRUE) {

  K <- get_attr(x, "K")
  p <- get_attr(x, "p")

  has_intercept <- get_attr(x, "intercept")
  add <- make_dummy_if(has_intercept)

  begin <- add + seq(1, p*K, K)
  end   <- add + seq(K, p*K, K)

  last <- nrow(x)
  if (last %ni% end)
    warning_glue("The partition did not divide the entire dataset")

  if (add == 1) {
    store <- vector("list", length = p + 1)
    store[[1]] <- x[1, , drop = FALSE]
    for (i in 1:p)
      store[[i + 1]] <- x[begin[i]:end[i],]
    num <- 0
  }else{
    store <- vector("list", length = p)
    for (i in 1:p)
      store[[i]] <- x[, begin[i]:end[i]]
    num <- 1
  }

  if (isTRUE(has_intercept) && isFALSE(keep_intercept))
    store <- store[-1]

  if (isFALSE(has_intercept) || isFALSE(keep_intercept)) {
    name_intercept <- NULL
  }else{
    name_intercept <- "intercept"
  }

  names(store) <- c(name_intercept, paste0(prefix_names, 1:p))

  if (isTRUE(transpose)) {
    purrr::map(store, t)
  }else{
    store
  }
}

# asfd --------------------------------------------------------------------



comp <- function(mat, arg_K = NULL, arg_p = NULL, intercept = NULL) {
  K <- arg_K %||% get_attr(mat, "K")
  p <- arg_p %||% get_attr(mat, "p")
  intercept <- intercept %||% get_attr(mat, "intercept")
  if (intercept) mat <- mat[-1,]
  pad <- eye(K * (p - 1), K*p )
  out <- rbind(t(mat), pad)
  meta_matrix(out, intercept = intercept, K = K, p = p)
}

uncomp <- function(mat, arg_K = NULL, arg_p = NULL) {
  K <- arg_K %||% get_attr(mat, "K")
  p <- arg_p %||% get_attr(mat, "p")
  t(mat[1:K,])
}


# Mimick some matlab functions --------------------------------------------

inv <- function(x) {
  solve(x)
}

assert_square_mat <- function(x) {
  # if (nrow(x) != ncol(x))
  #   stop("Non-square matrix")
}

vec <- function(x) {
  assert_square_mat(x)
  c(x)
}

vech <- function(x) {
  assert_square_mat(x)
  c(lower.tri(x, diag = TRUE))
}


eye <- function(n, p = NULL) {
  if (missing(p)) {
    diag(1, nrow = n, ncol = n)
  }else{
    diag(1, nrow = n, ncol = p)
  }
}

ones <- function(n, p = NULL) {
  if (missing(p)) {
    matrix(1, nrow = n, ncol = n)
  }else{
    matrix(1, nrow = n, ncol = p)
  }
}

zeros <- function(n, p = NULL) {
  if (missing(p)) {
    matrix(0, nrow = n, ncol = n)
  }else{
    matrix(0, nrow = n, ncol = p)
  }

}
