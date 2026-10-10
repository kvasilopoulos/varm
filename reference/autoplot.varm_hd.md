# Plot historical decompositions

Stacked columns of the shock contributions (and the baseline) to each
variable, with the actual series as a line.

## Usage

``` r
# S3 method for class 'varm_hd'
autoplot(object, responses = NULL, ...)

# S3 method for class 'varm_hd'
plot(x, ...)
```

## Arguments

- object, x:

  a `varm_hd` object from
  [`hd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md).

- responses:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variables to show.

- ...:

  unused.

## Value

A `ggplot` object.

## See also

Other plotting:
[`autoplot.varm_draws()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_draws.md),
[`autoplot.varm_forecast()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_forecast.md),
[`autoplot.varm_irf()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.md),
[`theme_varm()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)

## Examples

``` r
ggplot2::autoplot(hd(varm(yraw, p = 2)), responses = rate)
```
