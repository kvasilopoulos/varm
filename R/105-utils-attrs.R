

# Inherit attrs -----------------------------------------------------------

# TODO Crate matrix class that retain attributes
# meta_matrix

# Drop all previous attributes
set_attrs <- function(x, ...) {
  attrs <- dots_list(...)
  attributes(x) <- attrs
  x
}

inherit_attrs <- function(x, y) {
  attributes(x) <- attributes(y)
  x
}

add_attr <- function(x,  ...) {
  attrs <- dots_list(...)
  attributes(x) <- c(attributes(x), attrs)
  x
}

add_class <- function(x, nm = x) {
  class(x) <- append(class(x), nm)
  x
}

get_attr <- function(x, nm) {
  attr(x, nm)
}

attribute_names <- function(x) {
  names(attributes(x))
}

strip_attributes <- function(x, strip = NULL) {
  if (is.null(strip)) {
    exceptions <- c("names", "dimnames", "dim")
    idx_exceptions <- which(attribute_names(x) %in% exceptions )
    attributes(x)[-idx_exceptions] <- NULL
  }else{
    attributes(x)[strip] <- NULL
  }
  x
}

# mut_attr <- function(x, ...) {
#
#   nm_attr <- attr(x, nm)
#   if(is.null(nm_attr))
#     abort("Cannot mutate an non-existent attribute")
#
#   attr(x, nm) <- new_val
#   x
# }


