

# crossing ----------------------------------------------------------------
#similar to expand.grid

crossing(
  horizon = 1:10,
  impulse = c("infl", "ff"),
  response = c("infl", "ff")
)


# Analysis ----------------------------------------------------------------



library(tidyverse)
dt <- econdata::sw2001[,-1]

dt %>%
  abvar::spec(.endo_lags = 2) %>%
  abvar::var_ls() %>%
  summary()


object <- abvar::spec(dt, .endo_lags = 4) %>%
  var_ls() %>%
  irf(id = id_chol(), h =24, nboot = 100)

tidy.irf_var_ls(object) %>%
  ggplot() +
  geom_errorbar(aes(x = horizon, ymin = low, ymax = high)) +
  geom_point(aes(x = horizon, y = irf)) +
  facet_wrap(response ~ impulse, scales = "free_y") +
  theme_bw()


autoplot(object, primary_bands = bands(linetype = "dotted", opacity = 0.5)) +
  theme(panel.grid = element_blank())

tidy(object) %>%
  ggplot() +
  geom_line(aes(horizon, irf), col = "blue") +
  geom_line(aes(horizon, low), col = "red") +
  geom_line(aes(horizon, high), col = "red") +
  facet_wrap(~  impulse + response, scales = "free")

object <- obj <-  var_ls(abvar::spec(dt, .endo_lags = 4))

x1 <- irf(object, nboot = 200, h = 25, id = id_none())
autoplot(x1)

x2 <- irf(object, nboot = 200, h = 25, id = id_chol())
x3 <- irf(object, nboot = 200, h = 25, id = id_triangular())

autoplot(x1, primary_bands = null_bands()) +
  autolayer(x2, color = "blue")


irfs <- irf(object)

varnames <- attr(object, "endo_varnames")

filter_irf <- function(nms, x) {
  flt <- x %||% nms
}


impl <- varnames[which(varnames %in% impulse)]
resp <- varnames[which(varnames %in% response)]
tidy(irfs) %>%
  # filter(impulse %in% )
  mutate(title_names = paste(impulse, '\u2192', response)) %>%
  ggplot(aes(horizon, irf)) +
  geom_line() +
  facet_wrap(~title_names, scales = "free", dir = "h") +
  theme_default() +
  theme(
    strip.text.x = element_text(size = rel(1.2)),
    axis.title = element_blank()
  )


autoplot(x2, labeller = labeller(.multi_line = FALSE),  dir = "v",
         primary_bands = null_bands()) +
  autolayer(x1, color = "blue") +
  autolayer(x3)

autoplot(
  irfs, scales = "free", color = "blue", strip.position = "right",
  primary_bands = irf_bands(color = null_color(), shade_color = "#008080"),
  secondary_bands = irf_bands(probs = c(0.1, 0.9), color = null_color(), shade_color = "blue")
  ) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside")

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

autoplot(irfs, labeller = labeller(impulse = label_wrap_gen(10)),
         scales = "free", color = "blue",
  primary_bands = irf_bands(color = null_color(), shade_color = "#008080")) +


#
# cbind(
#   irfs$irfs[1,1,],
#   apply(irfs$boot_irfs, c(1,2,3), quantile, 0.16)[1,1,],
#   apply(irfs$boot_irfs, c(1,2,3), quantile, 0.84)[1,1,]
# ) %>%
#   plot.ts(plot.type = "single")
#
#
#



autoplot(irfs)
#
# c(irfs$irfs)
# matrix(irfs$irfs, ncol = 1)
#
# tidy(irfs) %>%
#   ggplot() +
#   geom_line(aes(horizon, irf)) +
#   geom_line(aes(horizon, low,)) +
#   geom_line(aes(horizon, high)) +
#   facet_wrap(impulse ~ response, scales = "free")
#

irfx <- irfs$irfs
array_to_tbl(irfx) %>%
  ggplot(aes(horizon, irf)) +
  geom_line() +
  facet_wrap(impulse ~ response, scales = "free_y")
autoplot(irfs)




econdata::sw2001[,-1] %>%
  abvar::spec(.endo_lags = 4) %>%
  var_ls() %>%
  irf(scale_factor= 1) %>%
  autoplot()


(x <-
  econdata::sw2001 %>%
  spec(infl + un + ff ~ 1) %>%
  var_ls() %>%
  irf())

econdata::sw2001 %>%
  spec(infl + un + ff ~ 1) %>%
  var_ols() %>%
  irf() %>%
  autoplot()

nboot <- 100

irf_boot <- function(x, nboot = 100) {

  A <- coefficients(x)
  p <- x$p
  K <- x$K
  h <- 12
  nboot = 100
  lhs <- as.matrix(x$lhs)

  if (is_bootstrapping){}
  resids <- residuals(x)
  resids_centered <- scale(resids, scale = FALSE)

  boot_irf <- array(0, dim = c(K, K, h + 1, nboot))
  for(i in 1:nboot) {

    resids_boot <- resids_centered[boot_inst(resids_centered), ]
    art_data <- t(A[-1,]) %>% tsDyn::VAR.boot(lag = 4, include = "none")
    # y_boot <- matrix(0, nrow = x$N, ncol = x$K)
    # lag_interval <- as.vector(as.matrix(x$lhs[1:p, ]))
    # y_boot[1:p, ] <- lag_interval
    #
    # for (t in 1:(x$N - x$p)) {
    #   y_boot[t + p, ] <- inv(A) %*% c(1, lag_interval) + resids_boot[t, ]
    #   add_row <- y_boot[t + p, ]
    # }
    boot_irf[,,,i] <- spec(as.data.frame(art_data)) %>% var_ls() %>% irf()
  }
}



# monte carlo -------------------------------------------------------------

boot_irf <- array(0, dim = c(K, K, h + 1, 1000))
for(i in 1:1000) {
  art_data <- t(A[-1,]) %>% VAR.sim(lag = 4, include = "none")
  boot_irf[,,,i] <- spec(as.data.frame(art_data)) %>% var_ls() %>% irf()
}


var_num <- 2
mpla <- data.frame(
  middle = irf(x)[var_num,,] %>% t(),
  lower = apply(boot_irf, c(1,2,3), quantile, probs = 0.05)[var_num,,] %>% t(),
  upper = apply(boot_irf, c(1,2,3), quantile, probs = 0.95)[var_num,,] %>% t(),
  horizon = 1:13
)


library(ggplot2)
ggplot(data = mpla, aes(x = horizon)) +
  geom_line(aes(y = middle.1)) +
  geom_line(aes(y = lower.1), col = "red") +
  geom_line(aes(y = upper.1), col = "red") +
  theme_bw()

