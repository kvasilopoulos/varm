# Tidy impulse responses and variance decompositions

Tidy impulse responses and variance decompositions

## Usage

``` r
# S3 method for class 'varm_irf'
tidy(x, bands = TRUE, ...)

# S3 method for class 'varm_fevd'
tidy(x, bands = TRUE, ...)

# S3 method for class 'varm_dyn'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  a `varm_irf` or `varm_fevd` object.

- bands:

  include the percentile bands (one row per cell and level).

- ...:

  unused.

- row.names, optional:

  unused.

## Value

a tibble with `response`, `shock` (factors in variable order), `horizon`
and `estimate`; with bands, additionally `level`, `conf.low`,
`conf.high` (one row per cell and level).

## See also

Other tidiers:
[`tidy.varm()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md),
[`tidy.varm_forecast()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_forecast.md),
[`tidy.varm_hd()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_hd.md)

## Examples

``` r
r <- irf(varm(yraw, p = 2), horizon = 4, ci = "mc", nboot = 50)
tidy(r)
#> # A tibble: 90 × 7
#>    response shock horizon estimate level conf.low conf.high
#>    <fct>    <fct>   <int>    <dbl> <dbl>    <dbl>     <dbl>
#>  1 infl     infl        0   0.308    0.9  0.279      0.334 
#>  2 unemp    infl        0  -0.0117   0.9 -0.0420     0.0285
#>  3 rate     infl        0   0.0877   0.9  0.00199    0.173 
#>  4 infl     unemp       0   0        0.9  0          0     
#>  5 unemp    unemp       0   0.332    0.9  0.310      0.359 
#>  6 rate     unemp       0  -0.268    0.9 -0.351     -0.197 
#>  7 infl     rate        0   0        0.9  0          0     
#>  8 unemp    rate        0   0        0.9  0          0     
#>  9 rate     rate        0   0.715    0.9  0.672      0.761 
#> 10 infl     infl        1   0.472    0.9  0.420      0.529 
#> # ℹ 80 more rows
tidy(r, bands = FALSE)
#> # A tibble: 45 × 4
#>    response shock horizon estimate
#>    <fct>    <fct>   <int>    <dbl>
#>  1 infl     infl        0   0.308 
#>  2 unemp    infl        0  -0.0117
#>  3 rate     infl        0   0.0877
#>  4 infl     unemp       0   0     
#>  5 unemp    unemp       0   0.332 
#>  6 rate     unemp       0  -0.268 
#>  7 infl     rate        0   0     
#>  8 unemp    rate        0   0     
#>  9 rate     rate        0   0.715 
#> 10 infl     infl        1   0.472 
#> # ℹ 35 more rows
```
