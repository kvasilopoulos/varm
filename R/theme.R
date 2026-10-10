#' varm plot theme and palette
#'
#' `theme_varm()` is the single theme used by every [ggplot2::autoplot()] method
#' in \pkg{varm}: black-and-white panel, dashed major grid, no strip background, legend at the bottom. `varm_palette()`
#' returns colours from the colour-blind safe Okabe-Ito palette.
#'
#' @param base_size base font size, passed to [ggplot2::theme_bw()].
#' @param n number of colours (at most 8; more are recycled with a warning).
#' @return `theme_varm()`: a [ggplot2::theme()]; `varm_palette()`: a character
#'   vector of hex colours.
#' @family plotting
#' @examples
#' varm_palette(3)
#' ggplot2::ggplot(yraw, ggplot2::aes(infl, rate)) +
#'   ggplot2::geom_point(colour = varm_palette(1)) +
#'   theme_varm()
#' @export
theme_varm <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linetype = "dashed")
    )
}

OKABE_ITO <- c(
  "#0072B2",
  "#D55E00",
  "#009E73",
  "#CC79A7",
  "#56B4E9",
  "#E69F00",
  "#F0E442",
  "#000000"
)

#' @rdname theme_varm
#' @export
varm_palette <- function(n = 8L) {
  if (n > length(OKABE_ITO)) {
    rlang::warn("`varm_palette()` has 8 colours; recycling")
  }
  rep_len(OKABE_ITO, n)
}
