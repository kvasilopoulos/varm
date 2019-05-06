
# Printing helpers --------------------------------------------------------

grey <- crayon::make_style("#949494")

gold <- crayon::make_style("#e6be8a")

# asfd --------------------------------------------------------------------


inv <- function(x) {
  solve(x)
}

# asdf --------------------------------------------------------------------


is_character0 <- function(x) {
  is.character(x) && length(x) == 0
}

`%ni%` <- Negate(`%in%`)


# defensive ---------------------------------------------------------------

warning_glue <- function(..., .sep = "", .envir = parent.frame(),
                         call. = FALSE, .domain = NULL) {
  warning(
    glue(..., .sep = .sep, .envir = .envir),
    call. = call., domain = .domain
  )
}

stop_glue <- function(..., .sep = "", .envir = parent.frame(),
                      call. = FALSE, .domain = NULL) {
  stop(
    glue(..., .sep = .sep, .envir = .envir),
    call. = call., domain = .domain
  )
}

