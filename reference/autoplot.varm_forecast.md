# Plot forecasts: history and fan chart

The last `n_history` observations are drawn as a line followed by the
mean forecast and nested interval ribbons.

## Usage

``` r
# S3 method for class 'varm_forecast'
autoplot(object, variables = NULL, n_history = 40L, ...)

# S3 method for class 'varm_forecast'
plot(x, ...)
```

## Arguments

- object, x:

  a `varm_forecast` object.

- variables:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variables to show.

- n_history:

  number of past observations to show.

- ...:

  unused.

## Value

A `ggplot` object.

## See also

Other plotting:
[`autoplot.varm_draws()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_draws.md),
[`autoplot.varm_hd()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_hd.md),
[`autoplot.varm_irf()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.md),
[`theme_varm()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)

## Examples

``` r
fc <- forecast(varm(yraw, p = 2), h = 8, nboot = 200)
ggplot2::autoplot(fc)

ggplot2::autoplot(fc, variables = rate, n_history = 20)
```
