# Plot local-projection results

[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
returns a `ggplot`;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints it and
returns it invisibly. Shaded ribbons are the stored band levels (widest
lightest).

## Usage

``` r
# S3 method for class 'lpm'
autoplot(object, ...)

# S3 method for class 'lp_state'
autoplot(object, type = c("states", "diff"), ...)

# S3 method for class 'lp_smooth'
autoplot(object, ...)

# S3 method for class 'lp_supt'
autoplot(object, ...)

# S3 method for class 'lp_forecast'
autoplot(object, ...)

# S3 method for class 'lpm'
plot(x, ...)

# S3 method for class 'lp_state'
plot(x, ...)

# S3 method for class 'lp_supt'
plot(x, ...)

# S3 method for class 'lp_forecast'
plot(x, ...)
```

## Arguments

- object, x:

  a result of
  [`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
  [`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
  [`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
  [`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
  [`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
  [`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)
  or
  [`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md).

- ...:

  unused ([`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  passes them to
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)).

- type:

  `lp_state` only: `"states"` or `"diff"`.

## Value

a `ggplot` object.

## Details

- `lpm`, `lp_joint`, `lp_panel`: response paths, one panel per response.

- `lp_state`: the high and low state overlaid (`type = "states"`) or
  their difference (`type = "diff"`).

- `lp_smooth`: smoothed against raw local projection.

- `lp_supt`: pointwise against sup-t simultaneous bands.

- `lp_forecast`: direct forecast fan chart.

## See also

Other local projections:
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
fit <- lp(yraw, shock = rate, horizon = 8, lags = 2, contemp = c(infl, unemp))
autoplot(fit)
```
