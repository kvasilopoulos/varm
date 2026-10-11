# Tidy and glance at local-projection results

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns one
row per horizon, response and (for
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md))
state, repeated for every confidence level so that bands are tidy too:
`horizon`, `response`, (`state`), `estimate`, `std.error`, `level`,
`conf.low`, `conf.high`.
[`glance()`](https://generics.r-lib.org/reference/glance.html) returns a
one-row summary of the fit. `as_tibble()` is an alias of
[`tidy()`](https://generics.r-lib.org/reference/tidy.html).

## Usage

``` r
# S3 method for class 'lpm'
tidy(x, conf.int = TRUE, ...)

# S3 method for class 'lp_state'
tidy(x, conf.int = TRUE, ...)

# S3 method for class 'lp_smooth'
tidy(x, conf.int = TRUE, raw = FALSE, ...)

# S3 method for class 'lp_forecast'
tidy(x, conf.int = TRUE, ...)

# S3 method for class 'lp_supt'
tidy(x, ...)

# S3 method for class 'lpm'
as_tibble(x, ...)

# S3 method for class 'lp_state'
as_tibble(x, ...)

# S3 method for class 'lp_forecast'
as_tibble(x, ...)

# S3 method for class 'lpm'
glance(x, ...)

# S3 method for class 'lp_panel'
glance(x, ...)

# S3 method for class 'lp_smooth'
glance(x, ...)

# S3 method for class 'lp_state'
glance(x, ...)
```

## Arguments

- x:

  a result of
  [`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
  [`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
  [`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
  [`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
  [`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
  [`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)
  or
  [`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md).

- conf.int:

  include the `level`, `conf.low` and `conf.high` columns (one row per
  band level).

- ...:

  unused.

- raw:

  `lp_smooth` only: also return the unsmoothed local projection as rows
  with `type == "raw"` (no bands).

## Value

a tibble.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)

## Examples

``` r
fit <- lp(yraw, shock = rate, horizon = 6, lags = 2, contemp = c(infl, unemp))
tidy(fit)
#> # A tibble: 42 × 7
#>    horizon response  estimate std.error level  conf.low conf.high
#>      <int> <fct>        <dbl>     <dbl> <dbl>     <dbl>     <dbl>
#>  1       0 infl     -1.24e-14  1.37e-13   0.9 -2.39e-13  2.14e-13
#>  2       1 infl     -5.89e- 2  4.20e- 2   0.9 -1.28e- 1  1.01e- 2
#>  3       2 infl     -6.16e- 2  7.23e- 2   0.9 -1.81e- 1  5.74e- 2
#>  4       3 infl     -3.15e- 2  1.05e- 1   0.9 -2.04e- 1  1.41e- 1
#>  5       4 infl     -2.17e- 2  1.52e- 1   0.9 -2.72e- 1  2.28e- 1
#>  6       5 infl     -7.74e- 3  1.71e- 1   0.9 -2.88e- 1  2.73e- 1
#>  7       6 infl     -1.53e- 2  1.79e- 1   0.9 -3.10e- 1  2.79e- 1
#>  8       0 unemp    -2.12e-14  1.23e-14   0.9 -4.14e-14 -9.19e-16
#>  9       1 unemp    -2.95e- 2  4.49e- 2   0.9 -1.03e- 1  4.43e- 2
#> 10       2 unemp    -7.90e- 2  7.56e- 2   0.9 -2.03e- 1  4.53e- 2
#> # ℹ 32 more rows
glance(fit)
#> # A tibble: 1 × 6
#>   shock horizon  lags se_type iv     nobs
#>   <chr>   <int> <int> <chr>   <lgl> <dbl>
#> 1 rate        6     2 nw      FALSE   213
```
