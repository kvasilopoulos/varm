#' Structural identification schemes
#'
#' Each constructor returns a small specification object consumed by [dynamics],
#' and friends. The reduced-form innovations are \eqn{u_t = P \varepsilon_t},
#' \eqn{E[\varepsilon_t\varepsilon_t']=I}; an identification scheme picks \eqn{P}
#' with \eqn{PP' = \Sigma}.
#'
#' * `id_chol()` recursive (Cholesky) ordering, Sims (1980).
#' * `id_longrun()` Blanchard-Quah (1989): the cumulative long-run effect matrix is lower triangular.
#' * `id_sign()` sign restrictions on impulse responses (Uhlig 2005; Rubio-Ramirez, Waggoner & Zha 2010).
#' * `id_given()` user supplied impact matrix \eqn{P}.
#'
#' @param order <[`tidy-select`][tidyselect::language]> variable ordering (first
#'   = most exogenous), e.g. `order = c(rate, infl, unemp)`; names, strings and
#'   positions all work. Defaults to the column order of the data.
#' @param restrictions a `K x K` matrix (rows: variables, columns: shocks) or a
#'   `K x K x H` array with entries `+1`, `-1` (sign) or `0` (unrestricted).
#'   Dimnames are used to label the shocks. Columns of all zeros are unidentified shocks.
#' @param horizon if `restrictions` is a matrix, the number of horizons (from impact)
#'   over which the restriction is imposed.
#' @param nrot number of accepted rotations to keep per reduced-form draw.
#' @param maxtry maximum rotations tried per reduced-form draw.
#' @param P a `K x K` impact matrix.
#' @family identification
#' @examples
#' fit <- varm(yraw, p = 2)
#' id_chol(order = c(rate, infl, unemp))
#' irf(fit, id_chol(order = c(rate, infl, unemp)), horizon = 4, ci = "none")
#'
#' S <- matrix(c(1, 1, 1, 1, -1, 0, 1, 0, 0), 3, 3, byrow = TRUE)
#' id_sign(S)
#' @name identification
NULL

new_id <- function(type, ..., class = paste0("id_", type)) {
  structure(list(type = type, ...), class = c(class, "varm_id"))
}

#' @rdname identification
#' @export
id_chol <- function(order = NULL) {
  new_id("chol", order = rlang::enquo(order))
}

#' @rdname identification
#' @export
id_longrun <- function() {
  new_id("longrun")
}

#' @rdname identification
#' @export
id_given <- function(P) {
  new_id("given", P = as.matrix(P))
}

#' @rdname identification
#' @export
id_sign <- function(restrictions, horizon = 1L, nrot = 1L, maxtry = 1000L) {
  R <- restrictions
  if (length(dim(R)) == 2L) {
    R <- array(
      rep(R, horizon),
      c(dim(R), horizon),
      dimnames = c(dimnames(R), list(NULL))
    )
  }
  if (length(dim(R)) != 3L || dim(R)[1] != dim(R)[2]) {
    rlang::abort("`restrictions` must be K x K or K x K x H")
  }
  R[is.na(R)] <- 0
  if (!all(R %in% c(-1, 0, 1))) {
    rlang::abort("restrictions must be -1, 0 or +1")
  }
  new_id("sign", S = R, nrot = as.integer(nrot), maxtry = as.integer(maxtry))
}

# printing ------------------------------------------------------------------

#' @export
print.varm_id <- function(x, digits = 3, ...) {
  switch(
    x$type,
    chol = print_id_chol(x),
    longrun = cat(
      "<varm_id> long-run (Blanchard-Quah)\n",
      "cumulative long-run effect matrix is lower triangular:\n",
      "shock j has no long-run effect on variable i for i < j\n",
      sep = ""
    ),
    given = {
      cat("<varm_id> user-supplied impact matrix P (response x shock)\n")
      print(round(x$P, digits))
    },
    sign = print_id_sign(x)
  )
  invisible(x)
}

print_id_chol <- function(x) {
  ord <- if (rlang::quo_is_null(x$order)) {
    "data column order"
  } else {
    rlang::as_label(x$order)
  }
  cat(
    "<varm_id> recursive (Cholesky), order: ",
    ord,
    "\n(first = most exogenous; the contemporaneous effect of later ",
    "variables on earlier ones is zero)\n",
    sep = ""
  )
}

# restriction matrix as +, -, and a dot for unrestricted
sign_symbols <- function(m) {
  out <- matrix("\u00b7", nrow(m), ncol(m))
  out[m > 0] <- "+"
  out[m < 0] <- "-"
  dimnames(out) <- list(
    rownames(m) %||% paste0("v", seq_len(nrow(m))),
    colnames(m) %||% paste0("e", seq_len(ncol(m)))
  )
  out
}

print_id_sign <- function(x) {
  S <- x$S
  H <- dim(S)[3]
  cat(
    "<varm_id> sign restrictions (",
    sum(S[,, 1] != 0),
    " at impact; ",
    H,
    " horizon",
    if (H > 1) "s",
    ")\n",
    "  + positive, - negative, \u00b7 unrestricted; rows = variables, ",
    "cols = shocks\n",
    "  rotations: ",
    x$nrot,
    " per draw, max ",
    x$maxtry,
    " tries\n",
    sep = ""
  )
  # group consecutive horizons with identical restrictions
  changed <- vapply(
    seq_len(H)[-1],
    function(h) !identical(S[,, h], S[,, h - 1L]),
    NA
  )
  grp <- cumsum(c(TRUE, changed))
  for (g in unique(grp)) {
    hs <- which(grp == g) - 1L
    if (max(grp) > 1L) {
      label <- if (length(hs) > 1) paste0(min(hs), "..", max(hs)) else hs
      cat("\nhorizon ", label, ":\n", sep = "")
    }
    print(sign_symbols(S[,, hs[1] + 1L]), quote = FALSE, right = TRUE)
  }
  free <- which(colSums(abs(S[,, 1, drop = FALSE])) == 0)
  if (length(free)) {
    shocks <- dimnames(S)[[2]] %||% paste0("e", seq_len(dim(S)[2]))
    cat("\nunrestricted shocks: ", paste(shocks[free], collapse = ", "), "\n")
  }
}

# Structural IRF draws [K,K,H+1,N'] and bookkeeping ---------------------------
structural_irf <- function(Bs, Ss, K, p, H, id, varnames) {
  if (is.function(id$structural)) {
    return(id$structural(Bs, Ss, K, p, H, id, varnames))
  }
  if (id$type == "sign") {
    return(sign_irf(Bs, Ss, K, p, H, id))
  }
  type <- c(chol = 0L, longrun = 1L, given = 2L)[[id$type]]
  ord <- chol_order(id, K, varnames, type)
  P <- if (type == 2L) id$P else matrix(0, K, K)
  cube <- cpp_map_irf(Bs, Ss, p, H, type, as.integer(ord - 1L), P)
  list(
    irf = as_irf4(cube, K, H + 1L),
    draw = seq_len(dim(Bs)[3]),
    accept_rate = NA_real_
  )
}

chol_order <- function(id, K, varnames, type) {
  ord <- select_vars(id$order, varnames, "variable in `order`")
  if (is.null(ord)) {
    return(seq_len(K))
  }
  if (type == 0L && (length(ord) != K || anyDuplicated(ord))) {
    rlang::abort("`order` must be a permutation of the variables")
  }
  ord
}

sign_irf <- function(Bs, Ss, K, p, H, id) {
  if (dim(id$S)[1] != K) {
    rlang::abort(
      "sign restrictions dimension does not match the number of variables"
    )
  }
  out <- cpp_sign_irf(Bs, Ss, p, H, id$S, id$nrot, id$maxtry)
  if (out$accepted == 0L) {
    rlang::abort(
      paste0(
        "no draw satisfied the sign restrictions (tries = ",
        out$tries,
        ")"
      )
    )
  }
  list(
    irf = as_irf4(out$irf, K, H + 1L),
    draw = out$draw,
    P = out$P,
    accept_rate = out$accepted / out$tries
  )
}
