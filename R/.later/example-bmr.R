library(BMR)
avar <- function(.data, lags = 4, irhor = 20) {
  var_data <- as.matrix(.data)
  var_obj <- new(cvar)
  var_obj$build(var_data, TRUE, lags)
  var_obj$estim()
  var_obj$boot(1000)
  plot_data <- irf_bmr(var_obj, irhor, percentiles = c(.16, .50, .84),
  )

  list(plot_data$plot_vals, colnames(.data))
}

irf_bmr <- function(obj, periods = 10, var_names = NULL,
                    percentiles = c(.05, .50, .95), which_shock = NULL, which_response = NULL) {
  if (periods <= 0) {
    stop("error: need periods to be > 0")
  }

  M <- obj$M
  n_draws <- dim(obj$beta_draws)[3]
  irf_temp <- obj$IRF(periods)$irf_vals
  # put the IRFs in a tesseract-type format

  irf_tess <- array(NA, dim = c(M, M, periods, n_draws))

  for (i in 1:n_draws) {
    irf_tess[, , , i] <- irf_temp[, , ((i - 1) * periods + 1):(i * periods)]
  }
  irf_tess_1 <- apply(irf_tess, c(3, 1, 2), sort) # fix annoying bug
  irf_tess <- aperm(irf_tess_1, c(2, 3, 1, 4))

  irf_upper <- min(round(percentiles[3] * n_draws), n_draws)
  irf_mid <- round(percentiles[2] * n_draws)
  irf_lower <- max(round(percentiles[1] * n_draws), 1)

  if (is.null(which_shock)) {
    which_shock <- 1:M
  }
  if (is.null(which_response)) {
    which_response <- 1:M
  }
  n_response <- length(which_response)
  n_shocks <- length(which_shock)

  if (class(var_names) != "character") {
    var_names <- character(length = M)
    for (i in 1:M) {
      var_names[i] <- paste("VAR", i, sep = "")
    }
  }

  plot_vals <- array(NA, dim = c(periods, 4, M, M))
  IRFPData <- 0
  for (i in 1:M) {
    for (k in 1:M) {
      IRFPData <- data.frame(irf_tess[, k, irf_lower, i], irf_tess[, k, irf_mid, i], irf_tess[, k, irf_upper, i], 1:(periods))
      IRFPData <- as.matrix(IRFPData)
      plot_vals[, , k, i] <- IRFPData
    }
  }
  return <- list(plot_vals = plot_vals)
}

plot_avar <- function(x, names = NULL, which_shock = 1, ncol = 4, magnitude = -1) {
  nms <- x[[2]]
  if (!is.null(names)) {
    nms <- names
  }
  lower <- x[[1]][, 1, , which_shock]
  mean <- x[[1]][, 2, , which_shock]
  upper <- x[[1]][, 3, , which_shock]

  gg <- list(lower, mean, upper) %>%
    map(mat_to_tbl, nms) %>%
    reduce(full_join, by = c("horizon", "name")) %>%
    set_names(c("horizon", "name", "lower", "irf", "upper")) %>%
    mutate_at(vars(upper, lower, irf), ~ magnitude*.x / irf[which_shock]) %>%
    ggplot() +
    geom_line(aes(horizon, irf)) +
    geom_line(aes(horizon, lower), linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(horizon, ymax = upper, ymin = lower), alpha = 0.2, fill = "grey50") +
    geom_line(aes(horizon, upper), linetype = "dashed", color = "grey50") +
    facet_wrap(~name, scales = "free", ncol = ncol) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.grid.major = element_blank()
      panel.grid.major = element_line(linetype = "dashed")
    )
  return(gg)
}

mat_to_tbl <- function(x, nms) {
  x %>%
    as_tibble() %>%
    `colnames<-`(nms) %>%
    add_column(horizon = 1:nrow(.)) %>%
    pivot_longer(-horizon, names_ptypes = list(name = factor(levels = nms)))
}

plot_avar(avar(dt), which_shock = 2)
