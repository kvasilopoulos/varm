

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


#'@export
bands <- function(probs = c(0.16, 0.84),
                  color = null_color(), linetype = "dashed", size = 0.6,
                  shade_color = "grey50", opacity = 0.2) {
  function(object, selected_impulses, selected_responses) {
    irf_data <- tidy(object, probs = probs) %>%
      filter(impulse %in% selected_impulses, response %in% selected_impulses) %>%
      mutate(imp_resp = forcats::as_factor(paste(impulse, "->", response)))
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
  function(object, selected_impulses, selected_responses) {
    irf_data <- tidy(object, probs = probs, sec_probs = prob_dodge) %>%
      filter(impulse %in% selected_impulses, response %in% selected_impulses) %>%
      mutate(imp_resp = forcats::as_factor(paste(impulse, "->", response)))
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
