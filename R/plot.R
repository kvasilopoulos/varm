#' Plot impulse responses / variance decompositions
#'
#' Grid of panels: rows = responding variable, columns = shock. Shaded bands are the
#' stored percentile bands (outermost level lightest).
#' @param x a `varm_irf` or `varm_fevd` object.
#' @param shocks,responses subsets (names or positions).
#' @param ... unused.
#' @return A `ggplot` object.
#' @return A `ggplot` object.
#' @importFrom rlang .data
#' @export
plot.varm_dyn <- function(x, shocks = NULL, responses = NULL, ...) {
  vn <- x$varnames
  K <- x$K
  sh <- if (is.null(shocks)) seq_len(K) else resolve_var(shocks, vn, "shock")
  rs <- if (is.null(responses)) {
    seq_len(K)
  } else {
    resolve_var(responses, vn, "response")
  }
  df <- as.data.frame(x)
  df <- df[df$shock %in% vn[sh] & df$response %in% vn[rs], ]
  df$shock <- factor(df$shock, vn[sh])
  df$response <- factor(df$response, vn[rs])
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$h, .data$estimate))
  nb <- if (is.null(x$bands)) 0L else dim(x$bands)[4] / 2
  if (nb) {
    qn <- grep("^q", names(df), value = TRUE)
    for (b in seq_len(nb)) {
      lo <- qn[b]
      hi <- qn[length(qn) + 1 - b]
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]),
          fill = "steelblue",
          alpha = 0.25
        )
    }
  }
  if (inherits(x, "varm_irf")) {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey50")
  }
  p +
    ggplot2::geom_line(colour = "navy", linewidth = 0.8) +
    ggplot2::facet_grid(
      response ~ shock,
      scales = "free_y",
      labeller = ggplot2::label_both
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw()
}

#' @rdname plot.varm_dyn
#' @param object a `varm_dyn` object.
#' @exportS3Method ggplot2::autoplot
autoplot.varm_dyn <- function(object, ...) plot.varm_dyn(object, ...)

#' @export
print.varm_dyn <- function(x, ...) {
  what <- if (inherits(x, "varm_fevd")) "FEVD" else "IRF"
  cat(
    "<varm ",
    what,
    "> K = ",
    x$K,
    ", horizon 0..",
    x$H,
    ", id = ",
    x$id$type,
    if (!is.null(x$bands)) {
      paste0(", bands: ", paste0(x$level * 100, "%", collapse = "/"))
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  if (!is.na(x$accept_rate %||% NA)) {
    cat(
      "sign-restriction acceptance rate: ",
      format(x$accept_rate, digits = 3),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Tidy data frame of an IRF / FEVD object
#' @param x a `varm_dyn` object.
#' @param row.names,optional unused.
#' @param ... unused.
#' @export
as.data.frame.varm_dyn <- function(x, row.names = NULL, optional = FALSE, ...) {
  g <- expand.grid(
    response = x$varnames,
    shock = x$varnames,
    h = 0:x$H,
    stringsAsFactors = FALSE
  )
  df <- cbind(g, estimate = as.vector(x$estimate))
  if (!is.null(x$bands)) {
    nm <- dimnames(x$bands)[[4]]
    for (b in seq_along(nm)) {
      df[[paste0("q", trimws(nm[b]))]] <- as.vector(x$bands[,,, b])
    }
  }
  df
}

#' @export
print.varm_forecast <- function(x, ...) {
  cat("<varm forecast> h = ", x$h, "\n", sep = "")
  print(round(x$mean, 4))
  invisible(x)
}
