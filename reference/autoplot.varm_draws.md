# Plot parameter draws: trace and density

Plot parameter draws: trace and density

## Usage

``` r
# S3 method for class 'varm_draws'
autoplot(object, type = c("trace", "density"), terms = NULL, ...)

# S3 method for class 'varm_draws'
plot(x, ...)
```

## Arguments

- object, x:

  a `varm_draws` object from
  [`draws()`](https://kvasilopoulos.github.io/varm/reference/draws.md)
  (or a Bayesian fit's draws).

- type:

  `"trace"` or `"density"`.

- terms:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  coefficient (row) names such as `infl.l1`, shown for every equation.
  Default: each variable's own first lag.

- ...:

  unused.

## Value

A `ggplot` object.

## See also

Other plotting:
[`autoplot.varm_forecast()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_forecast.md),
[`autoplot.varm_hd()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_hd.md),
[`autoplot.varm_irf()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.md),
[`theme_varm()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)

## Examples

``` r
d <- draws(varm(yraw, p = 2), "mc", n = 200)
ggplot2::autoplot(d)

ggplot2::autoplot(d, type = "density", terms = c(rate.l1, const))
```
