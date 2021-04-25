


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


# class and attrs ---------------------------------------------------------

set_class <- function(x, nm) {
  class(x) <- nm
  x
}

add_class <- function(x, ...) {
  class(x) <- append(c(...), class(x))
  x
}

set_attrs <- function(x, ...) {
  attrs <- dots_list(...)
  attributes(x) <- attrs
  x
}

#' @importFrom rlang dots_list
add_attr <- function(x,  ...) {
  attrs <- dots_list(...)
  attributes(x) <- c(attributes(x), attrs)
  x
}

inherit_attrs <- function(x, y) {

  attr_x <- attributes(x) %>% names() %||% NA_character_
  attr_y <- attributes(y) %>% names() %||% NA_character_

  remove_x <- which(attr_x %in% attr_y)
  attributes(y)[remove_x] <- NULL # remove duplicates

  attributes(x) <- c(attributes(x), attributes(y))
  x
}
