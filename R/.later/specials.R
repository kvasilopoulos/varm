# form <- function(rows = NULL, cols = NULL) {
#   is_rows_vars <- is.null(rows) || is_quosures(rows)
#     if (!is_rows_vars) {
#       if (!is.null(cols)) {
#         stop("`rows` must be `NULL` or a `vars()` list if `cols` is a `vars()` list", call. = FALSE)
#       }
#       # facets_list <- as_facets_list(rows)
#     }
#
# }
#
# # rlang
# library(rlang)
# a <- quote(a)
# b <- quote(b)
# new_formula(a+b, NULL)
#
# new_formula(quote(a), NULL)
# as_quosure(a ~ a)
#
# xreg <- function(...) {
#   new_formula(
#     lhs = purrr::reduce(enexprs(...), function(.x, .y) call2("+", .x, .y)),
#     rhs = 1
#   )
# }
#
# xreg(a,b)
#
# # fable -------------------------------------------------------------------
#
#
#
# model_xreg <- function(...){
#   model_formula <- new_formula(
#     lhs = NULL,
#     rhs = reduce(enexprs(...), function(.x, .y) call2("+", .x, .y))
#   )
#   env <- map(enquos(...), get_env)
#   env[map_lgl(env, compose(is_empty, env_parents))] <- NULL
#   env <- if(!is_empty(env)) get_env(env[[1]]) else base_env()
#   out <- model.frame(model_formula, data = env, na.action = stats::na.pass)
# }
#
#
# # trend -------------------------------------------------------------------
#
# new_specials <- function (..., .required_specials = NULL, .xreg_specials = NULL)
# {
#   specials <- squash(list2(...))
#   if (is.null(specials$xreg)) {
#     specials$xreg <- function(...) abort(sprintf("Exogenous regressors are not supported for %s.",
#                                                  self$model))
#   }
#   structure(specials, required_specials = .required_specials,
#             xreg_specials = .xreg_specials, class = "fable_specials")
# }
#
# trend.tbl_ts <- function(x, knots = NULL, origin = NULL){
#   idx_num <- x[[expr_text(tsibble::index(x))]] %>% units_since
#   knots_num <- if(is.null(knots)){NULL} else {knots %>% units_since}
#   index_interval <- interval(x) %>% time_unit()
#   idx_num <- idx_num/index_interval
#   knots_num <- knots_num/index_interval
#   if(!is.null(origin)){
#     origin <- units_since(origin)/index_interval
#   }
#
#   trend(idx_num, knots_num, origin)
# }
#
#
# trend.numeric <- function(x, knots = NULL, origin = NULL){
#   if(!is.null(origin)){
#     origin <- origin - 1 # trend should count from 1
#     x <- x - origin
#     knots <- knots - origin
#   }
#   knots_exprs <- knots %>%
#     map(function(.x) pmax(0, x-.x)) %>%
#     set_names(map_chr(knots, function(.x) paste0("trend_",format(.x))))
#   tibble(trend = x,
#          !!!knots_exprs)
# }
#
# common_xregs <- list(
#   trend = function(knots = NULL, origin = NULL){
#     if(is.null(origin)){
#       if(is.null(self$origin)){
#         self$origin <- self$data[[expr_text(index(self$data))]][[1]]
#       }
#       origin <- self$origin
#     }
#     fable:::trend(self$data, knots, origin) %>% as.matrix
#   }
# )
