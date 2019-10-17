irf.var_ls <- function(object, horizon = 12, nboot = 200, id = id_none(), scale = 1, ...) {

  K <- object$K
  p <- object$p
  A <- comp(object$coefficients, arg_K = K, arg_p = p, intercept = attr(object, "intercept"))
  B <- object$vcov
  Bid <- t(chol(B)) # id
  Bscale <- Bid/base::diag(Bid) * scale # shock_scale

  irfs <- irf_algo1(A, Bscale, K = K, p = p, h = horizon)
  boot_irfs <- irf_boot_var_ls(obj = object, h = horizon, nb = nboot)
  structure(
    list(
      irfs = irfs,
      boot_irfs = boot_irfs
    ),
    class = append("irf_var_ls", class(irfs))
  )
}

#' Impulse Response Function Algorithm 1
#'
#' @import expm expm
irf_algo1 <- function(A, B, K, p, h) {
  nms <- rownames(B)
  J <- cbind(eye(K), zeros(K, K*(p - 1)))
  irf_comp <-  array(0, c(K * p, h + 1, K * p))# matrix(NA, K^2, h + 1)
  irf_comp[1:K, 1, 1:K] <- B # diag(K) # B #  for unit shock scale_shock(id_shock(B))
  for (i in 1:h) {
    irf_comp[, i + 1, ] <- c(J %*% (A %^% i) %*% t(J) %*% B)
  }
  irf_comp <- irf_comp[1:K, , 1:K]
  out <- aperm(irf_comp, c(1, 3, 2))
  dimnames(out) <- list(nms, nms, NULL)
  out
}


#' boots <- irf_boot_var_ls(obj)
#'
irf_boot_var_ls <- function(obj, h, nb) {
  K <- obj$K
  p <- obj$p
  has_intercept <- get_attr(obj, "intercept")
  bty <- array(0, c(K, K, h + 1, nb))
  for (i in 1:nb) {
     auxy <- boot_var(obj)
     esty <- var_ls(abvar::spec(as.data.frame(auxy), .endo_lags = 4))
     A <- comp(esty$coefficients, arg_K = K, arg_p = p, intercept = has_intercept)
     B <- esty$vcov
     Bid <- t(chol(B))
     Bscale <- Bid/base::diag(Bid)
     bty[,,,i] <- irf_algo1(A, Bscale, K = K, p = p, h = h)
  }
  nms <- names(obj$lhs)
  dimnames(bty) <- list(nms, nms, NULL, NULL)
  bty
}

# se <- function() {}
# se_boot <- function(engine , nboot = 100, ci = c(0.05, 0.95), seed = NULLs) {}


scale_shock <- function(type = c("se", "unit"), factor = 1) {
  type <- match.arg(type)
  if (type == "se") {
    out <- function(x) x * factor
  }else{
    out <- function(x) 1/base::diag(x) * factor
  }
  out
}

#' @importFrom generics tidy
tidy.irf_var_ls <- function(x) {
  irfs <- array_to_tbl(x$irfs)
  boot_dist <- array_to_list(x$boot_irfs, margin = c(1, 2, 3)) %>%
    enframe() %>% select(distr = value)
  low <- apply(x$boot_irfs, c(1,2,3), quantile, probs = 0.16) %>%
    array_to_tbl() %>% select(low = irf)
  high <- apply(x$boot_irfs, c(1,2,3), quantile, probs = 0.84) %>%
    array_to_tbl() %>% select(high = irf)
  bind_cols(irfs, low, high) %>%
    bind_cols(boot_dist)
}

autoplot.irf_var_ls <- function(object, ...) {
  tidy(object) %>%
    # mutate(irf_names = paste(impulse, " ~ ", response)) %>%
    ggplot() +
    geom_line(aes(horizon, low), linetype = "dashed") +
    geom_line(aes(horizon, irf)) +
    geom_line(aes(horizon, high), linetype = "dashed") +
    facet_wrap(response ~ impulse, scales = "free_y") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      axis.title = element_blank()
    )
}
