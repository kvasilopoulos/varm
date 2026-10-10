# Plot impulse responses and variance decompositions

Grid of panels (rows: responding variable, columns: shock) with the
point or median response and nested ribbons for each stored level,
widest lightest.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
ggplot and returns it invisibly.

## Usage

``` r
# S3 method for class 'varm_irf'
autoplot(
  object,
  shocks = NULL,
  responses = NULL,
  ci = c("band", "line", "none"),
  compare = NULL,
  ...
)

# S3 method for class 'varm_fevd'
autoplot(
  object,
  shocks = NULL,
  responses = NULL,
  ci = c("band", "line", "none"),
  compare = NULL,
  ...
)

# S3 method for class 'varm_dyn'
plot(x, ...)
```

## Arguments

- object, x:

  a `varm_irf` or `varm_fevd` object.

- shocks, responses:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  subsets of the shocks and responding variables, e.g. `shocks = rate`.

- ci:

  `"band"` (shaded ribbons, default), `"line"` (dashed bounds) or
  `"none"`.

- compare:

  named list of further `varm_irf` / `varm_fevd` objects to overlay on
  `object` (labelled "Baseline" and the list names), one colour each.

- ...:

  unused.

## Value

A `ggplot` object.

## See also

Other plotting:
[`autoplot.varm_draws()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_draws.md),
[`autoplot.varm_forecast()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_forecast.md),
[`autoplot.varm_hd()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_hd.md),
[`theme_varm()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)

## Examples

``` r
fit <- varm(yraw, p = 2)
r <- irf(fit, horizon = 8, ci = "mc", nboot = 100)
ggplot2::autoplot(r)

ggplot2::autoplot(r, shocks = rate, responses = c(infl, unemp))

ggplot2::autoplot(r, ci = "line")

r2 <- irf(varm(yraw, p = 4), horizon = 8, ci = "none")
ggplot2::autoplot(r, compare = list("p = 4" = r2), ci = "none")
#> Error in rbind(deparse.level, ...): numbers of columns of arguments do not match
ggplot2::autoplot(fevd(fit, horizon = 8, ci = "none"))
```
