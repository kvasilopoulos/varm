
# Model_response ----------------------------------------------------------

#' @importFrom formula.tools lhs.vars
model_response <- function(.data, .formula) {
  # x <- f_lhs(fm)
  # ss <- character()
  # for(i in 2:(length(x))) ss[i-1] <- deparse(x[[i]])
  ss <- formula.tools::lhs.vars(.formula)
  .data[, ss]
}

# lag_matrix --------------------------------------------------------------

#' @importFrom purrr set_names map_df reduce
#' @importFrom dplyr bind_cols
#' @importFrom glue glue
tbl_lag <- function(x, dimension = 1) {
  temp <- vector("list", length = dimension)
  for (dims in 1:dimension) {
    temp[[dims]] <- map_df(x, lag, dims) %>%
      set_names(glue("{names(x)}_lag{dims}"))
  }
  reduce(temp, bind_cols)
}


mat_lag <- function(x, dimension = 1) {
  x %>%
    tbl_lag(dimension = 1) %>%
    slice(-c(1:dimension)) %>%
    as.matrix()
}

# Model Specification -----------------------------------------------------

#' Model Specification
#'
#' @importFrom modelr model_matrix
#' @importFrom rlang is_formula
#' @importFrom dplyr bind_cols
#' @examples
#'
#' .data <-
#'  tibble(
#'    x = rnorm(100),
#'    y = rnorm(100),
#'    z = rnorm(100)
#'  )
#'
#' .data %>%
#'   add_trend() %>%
#'   model_spec(x + y ~ 1) %>%
#' @export
spec <- function(.data, .formula, .endo_lags = 2, .exo_lags = 2) {

  .call <- match.call()

  if (missing(.formula)) {
    nms <- names(.data)
    formula_expr <- paste(paste(nms, collapse = " + "), "~", 1)
    .formula <- as.formula(formula_expr)
    .call$.formula <- .formula
  }else{
    isformula <- rlang::is_formula(.formula)
    if (!isformula) {
      stop_glue("provide a valid formula")
    }
  }

  .terms <- terms.formula(.formula)
  has_intercept <- rlang::as_logical(attr(.terms, "intercept"))
  has_trend <- "trend" %in% attr(.terms, "term.labels")

  lhs <- model_response(.data, .formula)
  endo_varnames <- names(lhs)

  rhs_basic <- model_matrix(.data, .formula) %>%
    when(has_intercept ~ rename(., Intercept = `(Intercept)`),~.)

  lhs_lags <- lhs %>% tbl_lag(.endo_lags)

  rhs <- rhs_basic %>%
    bind_cols(lhs_lags)

  exo_varnames <- rhs_basic %>%
    when(has_intercept ~ select(., -Intercept), ~.) %>%
    when(has_trend ~ select(., -trend), ~.) %>%
    names()

  # Order last
  if (!is_character0(exo_varnames)) {

    rhs <- rhs %>%
      select(-exo_varnames, everything())

    if (.exo_lags != 0) {
      exo_lag_matrix <- rhs_basic %>%
        select(exo_varnames) %>%
        tbl_lag(.exo_lags)

      rhs <- rhs %>%
        bind_cols(exo_lag_matrix)
    }
  }

  structure(
    tibble::lst(
      lhs = lhs,
      rhs = rhs,
      lhs_lags = lhs_lags,
      call = .call,
      endo_varnames = endo_varnames,
      intercept = has_intercept,
      trend = has_trend,
      exo_varnames = exo_varnames,
      num_endo = length(endo_varnames),
      num_exo = length(exo_varnames),
      endo_lags = .endo_lags,
      exo_lags = .exo_lags
    ),
    class = "model_spec"
  )
}


print.model_spec <- function(x) {


  cat(grey("## Model Specification\n"))

  cli::cat_line()
  cat("Call:\n")
  print(x$call)

  cli::cat_line()
  cat(grey(glue("# LHS: {NROW(x$lhs)} x {NCOL(x$lhs)} ~ ")))
  cat(paste(names(x$lhs), collapse = ", "))

  cli::cat_line()
  cat(grey(glue("# RHS: {NROW(x$rhs)} x {NCOL(x$rhs)} ~ ")))

  names_rhs <- names(x$rhs)
  if (x$intercept) {
    cat(gold("Intercept, "))
    names_rhs <- str_replace(names_rhs, "Intercept", "")
  }

  if (x$trend) {
    cat(gold("trend, "))
    names_rhs <- str_replace(names_rhs, "trend", "")
  }

  cat(paste(stringi::stri_remove_empty(names_rhs ), collapse = ", "))


}
