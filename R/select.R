# Tidyselect helpers ---------------------------------------------------------

# Wrap a plain value (names, positions) as a quosure; pass quosures through.
as_quo <- function(x) {
  if (rlang::is_quosure(x)) x else rlang::new_quosure(x, rlang::empty_env())
}

#' Resolve a tidyselect specification to integer positions
#'
#' Internal helper behind every argument of \pkg{varm} that names variables
#' (`order`, `cumulative`, `shock`, `responses`, ...). Accepts bare names,
#' strings, integer positions, `c()`, `all_of()`, `starts_with()`, ... evaluated
#' against `names`.
#'
#' @param quo a quosure (from [rlang::enquo()]) or a plain value.
#' @param names character vector of variable names to select from.
#' @param arg label used in error messages.
#' @return integer positions into `names` (`NULL` when `quo` is `NULL`).
#' @noRd
select_vars <- function(quo, names, arg = "variable") {
  quo <- as_quo(quo)
  if (rlang::quo_is_null(quo)) {
    return(NULL)
  }
  pos <- rlang::set_names(seq_along(names), names)
  out <- tryCatch(
    tidyselect::eval_select(quo, pos, allow_rename = FALSE),
    vctrs_error_subscript_oob = function(e) {
      rlang::abort(
        paste0("unknown ", arg, ": ", paste(e$i, collapse = ", ")),
        call = NULL
      )
    }
  )
  unname(as.integer(out))
}

# Same, but exactly one variable is required.
select_one <- function(quo, names, arg = "variable") {
  i <- select_vars(quo, names, arg)
  if (length(i) != 1L) {
    rlang::abort(paste0("`", arg, "` must select exactly one variable"))
  }
  i
}
