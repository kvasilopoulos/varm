# methods -----------------------------------------------------------------

#' @export
stable <- function(x) {
  Acomp <- comp(x$coefficients)
  roots <- eigen(Acomp)$values
  mod_roots <- Mod(roots)
  structure(
    tibble(eigen = roots, modulus = mod_roots),
    class = "stability"
  )
}

#' @export
#' @importFrom ggplot2 ggplot geom_path geom_point geom_segment coord_equal labs
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous
autoplot.stability <- function(x) {
  center <- c(0,0)
  diameter <- 2
  r <- diameter / 2
  t <- seq(0, 2 * pi, length.out = 1e3)
  x_points <- center[1] + r * cos(t)
  y_points <- center[2] + r * sin(t)
  xy_circle <- tibble(x = x_points, y = y_points)
  xy_points <- tibble(Im = Im(x$eigen), Re = Re(x$eigen))
  xy_points_outside <- filter(xy_points, Im > 1 | Re > 1)

  ggplot() +
    geom_path(data = xy_circle, aes(x, y)) +
    geom_point(data = xy_points, aes(Im, Re), color = "blue") +
    geom_segment(aes(x = -1, xend = 1, y = 0 , yend = 0)) +
    geom_segment(aes(x = 0, xend = 0, y = -1 , yend = 1)) +
    geom_point(data = xy_points_outside, aes(Im, Re), color = "red") +
    coord_equal() +
    labs(x = "Real", y = "Imaginary", titel = "Inverse AR roots") +
    scale_x_continuous(breaks = c(-1,0,1)) +
    scale_y_continuous(breaks = c(-1,0,1)) +
    theme_classic()
}


is_stable <- function(x) {
  model <- stable(x)
  if (any(model$modulus < 1))
    FALSE
  else TRUE
}

check_stable <- function(x) {
  if (!is_stable(x)) {
    warning("model is not stable")
  }
}
