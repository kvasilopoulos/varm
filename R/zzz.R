.onLoad <- function (libnames, pkgname) {
  op <- options()
  op.transx <- list(
  )
  toset <- !(names(op.transx) %in% names(op))
  if (any(toset))
    options(op.transx[toset])
  invisible()
}

# Set Global Variables to avoid NOTES in cmdchecks
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("impulse", "response", "horizon", "low", "high", "sec_low", "sec_low", "y")
  )
}
