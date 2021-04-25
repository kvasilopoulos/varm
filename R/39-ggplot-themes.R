
# https://github.com/jrnold/ggthemes/blob/master/R

theme_varm <- function() {
  theme(

  )
}

theme_default <- function(base_size = 11, base_family = "", base_line_size = base_size/22,
                          base_rect_size = base_size/22) {
  theme_grey(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size) %+replace%
    theme(
      strip.background = element_blank(),
      axis.title = element_blank(),
      strip.text.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "grey20"),
      panel.grid.major = element_line(colour = "grey92", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(fill = "white", colour = NA),
    )
}

#' @importFrom ggplot2 %+replace% theme_grey theme element_blank element_text rel
#' element_rect element_line
theme_basic <- function(base_size = 11, base_family = "", base_line_size = base_size/22,
                        base_rect_size = base_size/22) {
  theme_grey(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size) %+replace%
    theme(
      strip.background = element_blank(),
      axis.title = element_blank(),
      strip.text.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", size = rel(1)),
      # panel.border = element_rect(fill = NA, colour = "grey20"),
      panel.grid = element_blank(),
      # panel.grid = element_line(colour = "grey92"),
      # panel.grid.minor = element_line(size = rel(0.5)),
      legend.key = element_rect(fill = "white", colour = NA),
    )
}
