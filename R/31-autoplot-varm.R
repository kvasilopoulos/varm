#' @importFrom generics tidy
#' @import dplyr
#' @import tibble
#' @export
tidy.irf_varm <- function(x, probs = c(1, 0.9), sec_probs = c(0.16, 0.84)) {
  irfs <- array_to_tbl(x$irfs)
  boot_dist <- array_to_list(x$boot_irfs, margin = c(1, 2, 3)) %>%
    enframe() %>% select(distr = value)
  low <- apply(x$boot_irfs, c(1,2,3), quantile, probs = probs[1]) %>%
    array_to_tbl() %>% select(low = irf)
  high <- apply(x$boot_irfs, c(1,2,3), quantile, probs = probs[2]) %>%
    array_to_tbl() %>% select(high = irf)

  if (!is.null(sec_probs)) {
    sec_low <- apply(x$boot_irfs, c(1,2,3), quantile, probs = sec_probs[1]) %>%
      array_to_tbl() %>% select(sec_low = irf)
    sec_high <- apply(x$boot_irfs, c(1,2,3), quantile, probs = sec_probs[2]) %>%
      array_to_tbl() %>% select(sec_high = irf)
    out <- bind_cols(irfs, low, sec_low, high, sec_high) %>%
      bind_cols(boot_dist)
  }else{
    out <- bind_cols(irfs, low, high) %>%
      bind_cols(boot_dist)
  }
  out
}

# impulse = c(1,2,3)
# impulse = c(infl, un, ff)
# impulse = c("infl", "un", "ff")
#
select_names <- function(x, nms) {
  if(is.numeric(x)) {
    snames <- nms[x]
  }
  if(is.character(x)) {
    snames <- intersect(x, nms)
  }
  if(is.null(x)) {
    snames <- nms
  }
  snames
}







#'@export
#'@examples
#'
#'
autoplot.irf_varm <- function(object, impulse = NULL, response = NULL,
                              color = "black", linetype = "solid",
                              arrow_labels = TRUE, show_zero = FALSE,
                              primary_bands = bands(), secondary_bands = sec_bands(), ...) {
  dots <- rlang::dots_list(...)

  dnames <- dimnames(object$irfs)
  selected_impulses <- select_names(impulse, dnames[[1]])
  selected_responses <- select_names(response, dnames[[2]])

  facet_formula <- if(arrow_labels) ~ imp_resp else (impulse ~ response)
  facet_scales <- if (is.null(dots$scales)) "free" else NULL

  gg <- tidy(object, probs = NA, sec_probs = NA) %>%
    filter(impulse %in% selected_impulses, response %in% selected_impulses) %>%
    mutate(imp_resp = forcats::as_factor(paste(impulse, "->", response))) %>%
    ggplot() +
    geom_line(aes(horizon, irf), color = color, linetype = linetype) +
    facet_wrap(facet_formula, scales = facet_scales, ...) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_default()
  if(show_zero) {
    gg <- gg +
      geom_hline(yintercept = 0, linetype = "dashed")
  }

  gg  +
    primary_bands(object, selected_impulses, selected_responses) +
    secondary_bands(object, selected_impulses, selected_responses)
}



#'@export
autolayer.irf_varm <- function(object, color = "red", linetype = "solid",
                                 primary_bands = null_bands(),
                                 secondary_bands = null_bands(), ...) {
  layer_data <- tidy(object) %>%
    mutate(imp_resp = forcats::as_factor(paste(impulse, "->", response))) %>%
    mutate(irf = 1.1*irf)
  list(
    geom_line(aes(horizon, irf),
              color = color, linetype = linetype, data = layer_data),
    primary_bands(object),
    secondary_bands(object)
  )
}
