# Tidy historical decomposition

Tidy historical decomposition

## Usage

``` r
# S3 method for class 'varm_hd'
tidy(x, ...)
```

## Arguments

- x:

  a `varm_hd` object from
  [`hd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md).

- ...:

  unused.

## Value

a long tibble with `time`, `response`, `shock`, `contribution` and
`actual`. The deterministic/initial-condition part appears as
`shock == "baseline"` rows (when available).

## See also

Other tidiers:
[`tidy.varm()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md),
[`tidy.varm_forecast()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_forecast.md),
[`tidy.varm_irf()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_irf.md)

## Examples

``` r
tidy(hd(varm(yraw, p = 2)))
#> # A tibble: 2,556 × 5
#>     time response shock contribution actual
#>    <int> <fct>    <fct>        <dbl>  <dbl>
#>  1     1 infl     infl        -0.498  1.01 
#>  2     2 infl     infl        -0.752  0.981
#>  3     3 infl     infl        -0.568  1.33 
#>  4     4 infl     infl        -0.329  1.30 
#>  5     5 infl     infl        -0.345  0.909
#>  6     6 infl     infl        -0.382  0.709
#>  7     7 infl     infl        -0.330  0.685
#>  8     8 infl     infl        -0.127  1.06 
#>  9     9 infl     infl         0.474  1.84 
#> 10    10 infl     infl         0.674  2.39 
#> # ℹ 2,546 more rows
```
