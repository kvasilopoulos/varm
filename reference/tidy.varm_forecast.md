# Tidy forecasts

Tidy forecasts

## Usage

``` r
# S3 method for class 'varm_forecast'
tidy(x, ...)
```

## Arguments

- x:

  a `varm_forecast` object from
  [`forecast()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md).

- ...:

  unused.

## Value

a tibble with `horizon`, `variable`, `estimate` (the mean forecast) and,
when intervals exist, `level`, `conf.low`, `conf.high` (one row per cell
and level).

## See also

Other tidiers:
[`tidy.varm()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md),
[`tidy.varm_hd()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_hd.md),
[`tidy.varm_irf()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_irf.md)

## Examples

``` r
tidy(forecast(varm(yraw, p = 2), h = 4, nboot = 100))
#> # A tibble: 24 × 6
#>    horizon variable estimate level conf.low conf.high
#>      <int> <fct>       <dbl> <dbl>    <dbl>     <dbl>
#>  1       1 infl         2.72   0.9     2.27      3.27
#>  2       2 infl         2.66   0.9     1.79      3.48
#>  3       3 infl         2.64   0.9     1.33      3.93
#>  4       4 infl         2.66   0.9     1.12      4.37
#>  5       1 unemp        4.92   0.9     4.33      5.41
#>  6       2 unemp        5.05   0.9     3.98      5.95
#>  7       3 unemp        5.16   0.9     3.99      6.35
#>  8       4 unemp        5.23   0.9     4.00      6.43
#>  9       1 rate         4.62   0.9     3.64      6.08
#> 10       2 rate         4.45   0.9     2.96      6.14
#> # ℹ 14 more rows
```
