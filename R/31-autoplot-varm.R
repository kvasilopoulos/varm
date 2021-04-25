#' @importFrom generics tidy
#' @importFrom dplyr full_join select
#' @export
tidy.irf_varm <- function(x, probs = c(0.1, 0.9), sec_probs = c(0.16, 0.84), ...) {


  # boot_dist <- array_to_list(x$boot_irfs, margin = c(1, 2, 3)) %>%
  #   enframe() %>% select(distr = value)
  # low <- apply(x$boot_irfs, c(1,2,3), quantile, probs = probs[1]) %>%
    # array_to_tbl() %>% select(low = irf)
  # high <- apply(x$boot_irfs, c(1,2,3), quantile, probs = probs[2]) %>%
  # array_to_tbl() %>% select(high = irf)

  irfs <- array_to_tbl(x$irfs)
  boot_dist <- array_to_list(x$boot_irfs, margin = c(1,2,3))

  full_join(irfs, boot_dist, by = c("impulse", "response", "horizon")) %>%
    # add_column(low) %>%
    update_tidy(probs = probs, sec_probs = sec_probs) %>%
    select(which(colMeans(!is.na(.)) == 1))
}

#' @importFrom dplyr mutate
#' @importFrom purrr map_dbl
update_tidy <- function(x, probs = c(0.1, 0.9), sec_probs = c(0.16, 0.84)) {
  x %>%
    mutate(
      low = map_dbl(distr, quantile, probs[1]),
      high = map_dbl(distr, quantile, probs[2]),
      sec_low = map_dbl(distr, quantile, sec_probs[1]),
      sec_high = map_dbl(distr, quantile, sec_probs[2])
    )
}

grid_names <- function(x) {
  dims <- dim(x)
  nms <- rownames(x)
  expand.grid(
    impulse = nms,
    response =  nms,
    horizon = 1:dims[3])
}

#' 3d array into tibble
#' @importFrom tidyr expand_grid
array_to_tbl <- function(x) {
  as_tibble(grid_names(x)) %>%
    add_column(irf = c(x)) %>% # c() to collapse dims of x
    arrange(impulse, response, horizon)
}

#' @importFrom purrr map
#' @importFrom dplyr arrange
#' @importFrom tibble as_tibble add_column
array_to_list <- function(x, margin = 3) {
  out <- apply(x, margin, list)
  dm <- dim(x)[4]
  boot <- map(out, unlist) %>%
    map(set_class, glue("boot"))
  as_tibble(grid_names(x)) %>%
    add_column(distr = boot) %>% # c() to collapse dims of x
    arrange(impulse, response, horizon)
}



# plotting ----------------------------------------------------------------



#' impulse = c(1,2,3)
#' impulse = c(infl, un, ff)
#'
select_names <- function(x, nms = var_names) {
  if(is.numeric(x)) {
    snames <- nms[x]
  }
  if(is.character(x)) {
    stopifnot(x %in% nms)
    snames <- nms
  }
  if(is.null(x)) {
    snames <- nms
  }
  snames
}

symbols <- list()
symbols$rarr <- "<span style='font-size:16px'> &rarr; </span>"

fcol <- function(x, color = "red") {
  glue("<span style='color:{color};'> {x} </span>")
}



# text <- c(
#   "Some text **in bold.**",
#   "Linebreaks<br>Linebreaks<br>Linebreaks",
#   "*x*<sup>2</sup> + 5*x* + *C*<sub>*i*</sub>",
#   "Some <span style='color:blue'>blue text **in bold.**</span><br>And *italics text.*<br>And some <span style='font-size:18pt; color:black'>large</span> text."
# )

#' Plotting methods for objects of `varm`
#'
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs theme geom_hline
#' @importFrom ggtext element_markdown
#' @importFrom dplyr filter mutate
#' @export
#' @examples
#'
#'
autoplot.irf_varm <- function(object, impulse = NULL, response = NULL, var_names = NULL,
                              color = "black", linetype = "solid", show_zero = FALSE,
                              label = "{fcol(impulse, 'red')} to {fcol(response, 'blue')}",
                              primary_bands = bands(), secondary_bands = sec_bands(), ...) {
  dots <- rlang::dots_list(...)

  dnames <- dimnames(object$irfs)
  selected_impulses <- select_names(impulse, var_names %||% dnames[[1]])
  selected_responses <- select_names(response, var_names %||% dnames[[2]])

  # TODO check label

  # facet_formula <- if(show_labels) ~ labels else (impulse ~ response)
  facet_scales <- if (is.null(dots$scales)) "free" else NULL

  tidy_object <- object %>%
    set_varnames(var_names) %>%
    tidy(probs = NA, sec_probs = NA) %>%
    filter(impulse %in% selected_impulses, response %in% selected_impulses) %>%
    mutate(labels = forcats::as_factor(glue::glue(label)))
  # return(tidy_object)
  gg <-  ggplot(tidy_object) +
    geom_line(aes(horizon, irf), color = color, linetype = linetype) +
    facet_wrap(~labels, scales = facet_scales) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_default() +
    theme(
      strip.text.x = element_markdown() # enable md in labels
    )
  if(show_zero) {
    gg <- gg +
      geom_hline(yintercept = 0, linetype = "dashed")
  }
  gg  +
    primary_bands(tidy_object) +
    secondary_bands(tidy_object)
}


split_plot <- function(x) {
  facets <- x$facet$params$facets
  facet_names <- names(facets)
  facet_vars <- unique(x$data[facet_names])[[facet_names]]
  nfacets <- NROW(facet_vars)
  free <- "free"
  map(1:nfacets, ~
        x + ggforce::facet_wrap_paginate(facets, scales = free, nrow = 1, ncol = 1, page = .x)
  ) %>%
    set_names(facet_vars)
}

# splitFacet <- function(x){
#   facet_vars <- names(x$facet$params$facets)         # 1
#   x$facet    <- ggplot2::ggplot()$facet              # 2
#   datasets   <- split(x$data, x$data[facet_vars])    # 3
#   new_plots  <- lapply(datasets,function(new_data) { # 4
#     x$data <- new_data
#     x})
# }


#https://stackoverflow.com/questions/30510898/split-facet-plot-into-list-of-plots
# splitFacet <- function(x, n = NULL){
#   facet_vars <- names(x$facet$params$facets)               # 1
#   if(is.null(n)){
#     x$facet  <- ggplot2::ggplot()$facet                    # 2a
#     datasets <- split(x$data, x$data[facet_vars])          # 3a
#   } else {
#     inter0 <- interaction(x$data[facet_vars], drop = TRUE) # 2b
#     inter  <- ceiling(as.numeric(inter0)/n)
#     datasets <- split(x$data, inter)                       # 3b
#   }
#   new_plots  <- lapply(datasets,function(new_data) {       # 4
#     x$data <- new_data
#     x})
# }



#' @importFrom ggplot2 geom_line
#'@export
autolayer.irf_varm <- function(object, color = "red", linetype = "solid",
                               primary_bands = null_bands(),
                               secondary_bands = null_bands(), ...) {
  layer_data <- tidy(object) %>%
    mutate(imp_resp = forcats::as_factor(paste(impulse, "->", response))) %>%
    mutate(irf = 1.1*irf)
  list(
    geom_line(aes(horizon, irf), color = color, linetype = linetype, data = layer_data),
    primary_bands(object),
    secondary_bands(object)
  )
}




#'@export
null_bands <- function() {
  function(object) {
    NULL
  }
}

#'@export
null_color <- function() {
  "#ffffff00"
}


#' @importFrom ggplot2 geom_line geom_ribbon
#'@export
bands <- function(probs = c(0.16, 0.84),
                  color = null_color(), linetype = "dashed", size = 0.6,
                  shade_color = "grey50", opacity = 0.2) {

  function(object) {

    irf_data <- update_tidy(object, probs = probs)

    list(
      geom_line(
        aes(horizon, low), data = irf_data,
        size = size, color = color, linetype = linetype),
      geom_line(
        aes(horizon, high), data = irf_data,
        size = size, color = color, linetype = linetype),
      geom_ribbon(
        aes(x = horizon, ymax = high, ymin = low), data = irf_data,
        fill = shade_color, alpha = opacity)
    )
  }
}

#'@export
sec_bands <- function(probs = c(0.1, 0.9), prob_dodge = c(0.16, 0.84),
                      color = null_color(), linetype = "dashed", size = 0.6,
                      shade_color = "grey75", opacity = 0.3) {
  function(object) {

    irf_data <- update_tidy(object, probs = probs, sec_probs = prob_dodge)

    list(
      geom_line(
        aes(horizon, low), data = irf_data,
        size = size, color = color, linetype = linetype),
      geom_line(
        aes(horizon, high), data = irf_data,
        size = size, color = color, linetype = linetype),
      geom_ribbon(
        aes(x = horizon, ymax = high, ymin = sec_high), data = irf_data,
        fill = shade_color, alpha = opacity),
      geom_ribbon(
        aes(x = horizon, ymax = sec_low, ymin = low), data = irf_data,
        fill = shade_color, alpha = opacity)
    )
  }
}
