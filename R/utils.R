
#' 3d array into tibble
#' @importFrom tidyr expand_grid
array_to_tbl <- function(x) {
  dims <- dim(x)
  nms <- rownames(x)
  grid_names <- expand.grid(
    impulse = nms,
    response =  nms,
    horizon = 1:dims[3])
  grid_names %>%
    add_column(irf = c(x)) %>% # c() to collapse dims of x
    arrange(impulse, response, horizon)
}

# custom colors -----------------------------------------------------------

grey <- crayon::make_style("#949494")

gold <- crayon::make_style("#e6be8a")


# helpers ----------------------------------------------------------------


is_character0 <- function(x) {
  is.character(x) && length(x) == 0
}

`%ni%` <- Negate(`%in%`)

`%||%` <- function(x, y) {
  if (is.null(x))
    y
  else x
}

meta_matrix <- function(x, ...) {
  structure(x, ...)
}

make_dummy_if <- function(condition) {
  if (condition) {
    1
  } else{
    0
  }
}



