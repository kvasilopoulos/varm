
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
#' @importFrom dplyr bind_cols lag
#' @importFrom glue glue
tbl_lag <- function(x, dimension = 1) {
  temp <- vector("list", length = dimension)
  for (dims in 1:dimension) {
    temp[[dims]] <- map_df(x, dplyr::lag, dims) %>%
      set_names(glue("{names(x)}_lag{dims}"))
  }
  reduce(temp, bind_cols)
}


mat_lag <- function(x, lags = 1) {
  dt <- as.matrix(x)
  nc <- ncol(dt)
  embed(dt, dimension = lags + 1)[, -1, drop = FALSE]
}


# trend -------------------------------------------------------------------

#' @importFrom rlang eval_bare enexpr
custom_trend <- function(n, .f) {
  if (missing(.f)) {
    .expr <- rlang::expr(.t)
  }else{
    .expr <- rlang::enexpr(.f)
  }
  .env <- rlang::env(.t = 1:n)
  rlang::eval_bare(.expr, env = .env)
}

trend_fun <- function(type = c("linear", "quadratic")) {
  function(x) {
    if (type == "linear") {
      fun <- custom_trend(x, .t)
    }else{
      fun <- custom_trend(x, .t + .t^2)
    }
    fun
  }
}


# Model Specification -----------------------------------------------------

# TODO parse the index
# TODO accomodate different specifications such as favar tvp
# set_prior
# TODO spec can be optional for more detailed specification and then provided for estimation
# alternatively avar


#' Create a model specification that can fit into multiple models
#'
#' @param .data df
#' @param .formula the formula
#' @param .endo_lags  number of endogenous lags
#' @param .exo_lags number of exogenous lags
#'
#' @importFrom modelr model_matrix
#' @importFrom rlang is_formula
#' @importFrom dplyr bind_cols everything select tibble slice rename
#' @importFrom purrr when
#' @importFrom stats as.formula terms.formula
#'
#' @export
#' @examples
#' spec(yraw)
#'
spec <- function(
  .data,
  .formula = NULL, # if null no formula is used
  .index = NULL, # TODO parse index
  .ordering = NULL, # TODO rearrangee for selection and/or ordering
  .endo_lags = 2,
  .endo_names = NULL,
  .trend = FALSE,
  .factor_names = NULL,
  .tvp = NULL,
  .sv = NULL,
  .exo_names = NULL,
  .exo_lags = NULL
) {
# TODO global spec to be used in every function and passed args through ellipsis
  nr <- nrow(.data)

  if (is.character(.trend)) {
    type <- match.arg(.trend, choices = c("linear", "quadratic"))
    trend_fun(type)(nr)
  }
  if (is.function(.trend)) {
  }

  .call <- match.call()

  # TODO handle index in spec
  # TODO sort out which arguments need a .

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
  has_trend <- if (is.null(.trend)) TRUE else FALSE
    #"trend" %in% attr(.terms, "term.labels")

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
    list(
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
    class = "spec"
  )
}

# ?pillar::dim_desc(dt)
# TODO add index parsing

#' @importFrom crayon cyan
#' @importFrom stringr str_replace
#' @export
print.spec <- function(x, ...) {
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
    cat(crayon::cyan("Intercept, "))
    names_rhs <- str_replace(names_rhs, "Intercept", "")
  }
  if (x$trend) {
    cat(gold("trend, "))
    names_rhs <- str_replace(names_rhs, "trend", "")
  }
  cat(paste(stringi::stri_remove_empty(names_rhs ), collapse = ", "))
}



# Specify prior -----------------------------------------------------------



prior_spec <- function(type) {

}

prior_analytical <- function() {

}

prior_minn <- function(coef = NULL, variance = NULL, decay = NULL, h1 = NULL,
                       h2 = NULL, h3 = NULL, h4 =NULL) {

}
