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
#' @param order variable ordering, as names or integer positions (first = most exogenous);
#'   defaults to the column order of the data.
#' @param restrictions a `K x K` matrix (rows: variables, columns: shocks) or a
#'   `K x K x H` array with entries `+1`, `-1` (sign) or `0` (unrestricted).
#'   Dimnames are used to label the shocks. Columns of all zeros are unidentified shocks.
#' @param horizon if `restrictions` is a matrix, the number of horizons (from impact)
#'   over which the restriction is imposed.
#' @param nrot number of accepted rotations to keep per reduced-form draw.
#' @param maxtry maximum rotations tried per reduced-form draw.
#' @param P a `K x K` impact matrix.
#' @name identification
NULL

#' @rdname identification
#' @export
id_chol <- function(order = NULL) {
  structure(list(type = "chol", order = order), class = c("id_chol", "varm_id"))
}

#' @rdname identification
#' @export
id_longrun <- function() {
  structure(list(type = "longrun"), class = c("id_longrun", "varm_id"))
}

#' @rdname identification
#' @export
id_given <- function(P) {
  structure(
    list(type = "given", P = as.matrix(P)),
    class = c("id_given", "varm_id")
  )
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
    stop("`restrictions` must be K x K or K x K x H", call. = FALSE)
  }
  R[is.na(R)] <- 0
  if (!all(R %in% c(-1, 0, 1))) {
    stop("restrictions must be -1, 0 or +1", call. = FALSE)
  }
  structure(
    list(
      type = "sign",
      S = R,
      nrot = as.integer(nrot),
      maxtry = as.integer(maxtry)
    ),
    class = c("id_sign", "varm_id")
  )
}

#' @export
print.varm_id <- function(x, digits = 3, ...) {
  grid <- function(m) {
    # restriction matrix, rows = variables, cols = shocks
    dn <- dimnames(m)
    rn <- dn[[1]] %||% paste0("v", seq_len(nrow(m)))
    cn <- dn[[2]] %||% paste0("e", seq_len(ncol(m)))
    print(
      matrix(m, nrow(m), dimnames = list(rn, cn)),
      quote = FALSE,
      right = TRUE
    )
  }
  switch(
    x$type,
    chol = cat(
      "<varm_id> recursive (Cholesky), order: ",
      if (is.null(x$order)) {
        "data column order"
      } else {
        paste(x$order, collapse = " -> ")
      },
      "\n(first = most exogenous; the contemporaneous effect of later variables on earlier ones is zero)\n",
      sep = ""
    ),
    longrun = cat(
      "<varm_id> long-run (Blanchard-Quah)\ncumulative long-run effect matrix is lower triangular:\n",
      "shock j has no long-run effect on variable i for i < j\n",
      sep = ""
    ),
    given = {
      cat("<varm_id> user-supplied impact matrix P (response x shock)\n")
      print(round(x$P, digits))
    },
    sign = {
      S <- x$S
      H <- dim(S)[3]
      sym <- function(m) {
        o <- matrix("\u00b7", nrow(m), ncol(m), dimnames = dimnames(m)[1:2])
        o[m > 0] <- "+"
        o[m < 0] <- "-"
        o
      }
      cat(
        "<varm_id> sign restrictions (",
        sum(S[,, 1] != 0),
        " at impact; ",
        H,
        " horizon",
        if (H > 1) "s",
        ")\n",
        "  + positive, - negative, \u00b7 unrestricted; rows = variables, cols = shocks\n",
        "  rotations: ",
        x$nrot,
        " per draw, max ",
        x$maxtry,
        " tries\n",
        sep = ""
      )
      # group consecutive horizons with identical restrictions
      new <- c(
        TRUE,
        vapply(
          seq_len(H)[-1],
          function(h) !identical(S[,, h], S[,, h - 1L]),
          NA
        )
      )
      grp <- cumsum(new)
      for (g in unique(grp)) {
        hs <- which(grp == g) - 1L
        if (max(grp) > 1L) {
          cat(
            "\nhorizon ",
            if (length(hs) > 1) paste0(min(hs), "..", max(hs)) else hs,
            ":\n",
            sep = ""
          )
        }
        grid(sym(S[,, hs[1] + 1L]))
      }
      free <- which(colSums(abs(S[,, 1, drop = FALSE])) == 0)
      if (length(free)) {
        cat(
          "\nunrestricted shocks: ",
          paste(
            (dimnames(S)[[2]] %||% paste0("e", seq_len(dim(S)[2])))[free],
            collapse = ", "
          ),
          "\n",
          sep = ""
        )
      }
    }
  )
  invisible(x)
}

# Structural IRF draws [K,K,H+1,N'] and bookkeeping ---------------------------
structural_irf <- function(Bs, Ss, K, p, H, id, varnames) {
  H1 <- H + 1L
  if (id$type == "sign") {
    if (dim(id$S)[1] != K) {
      stop(
        "sign restrictions dimension does not match the number of variables",
        call. = FALSE
      )
    }
    out <- cpp_sign_irf(Bs, Ss, p, H, id$S, id$nrot, id$maxtry)
    if (out$accepted == 0L) {
      stop(
        "no draw satisfied the sign restrictions (tries = ",
        out$tries,
        ")",
        call. = FALSE
      )
    }
    return(list(
      irf = as_irf4(out$irf, K, H1),
      draw = out$draw,
      P = out$P,
      accept_rate = out$accepted / out$tries
    ))
  }
  type <- c(chol = 0L, longrun = 1L, given = 2L)[[id$type]]
  ord <- if (is.null(id$order)) {
    seq_len(K)
  } else {
    resolve_var(id$order, varnames, "variable in `order`")
  }
  if (type == 0L && (length(ord) != K || anyDuplicated(ord))) {
    stop("`order` must be a permutation of the variables", call. = FALSE)
  }
  P <- if (type == 2L) id$P else matrix(0, K, K)
  cube <- cpp_map_irf(Bs, Ss, p, H, type, as.integer(ord - 1L), P)
  list(
    irf = as_irf4(cube, K, H1),
    draw = seq_len(dim(Bs)[3]),
    accept_rate = NA_real_
  )
}
