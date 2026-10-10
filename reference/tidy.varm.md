# Tidy a VAR fit

Coefficient table of a
[`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md) fit
in the broom convention, a one-row model summary, and the data augmented
with fitted values and residuals.

## Usage

``` r
# S3 method for class 'varm'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'varm'
glance(x, ...)

# S3 method for class 'varm'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  a `varm` fit.

- conf.int:

  include `conf.low` / `conf.high` (t-based).

- conf.level:

  confidence level for `conf.int`.

- ...:

  unused.

- data:

  data to augment; defaults to the data used for the fit.

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html): a tibble
with `response`, `term`, `estimate`, `std.error`, `statistic`, `p.value`
(and `conf.low`, `conf.high`).
[`glance()`](https://generics.r-lib.org/reference/glance.html): one row
with `K`, `p`, `nobs`, `logLik`, `AIC`, `BIC`, `HQ`, `max_root`.
[`augment()`](https://generics.r-lib.org/reference/augment.html): the
data plus `.fitted_<var>` and `.resid_<var>` columns (`NA` for the first
`p` rows).

## See also

Other tidiers:
[`tidy.varm_forecast()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_forecast.md),
[`tidy.varm_hd()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_hd.md),
[`tidy.varm_irf()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_irf.md)

## Examples

``` r
fit <- varm(yraw, p = 2)
tidy(fit)
#> # A tibble: 21 × 6
#>    response term     estimate std.error statistic  p.value
#>    <chr>    <chr>       <dbl>     <dbl>     <dbl>    <dbl>
#>  1 infl     infl.l1   1.54       0.0596    25.8   2.13e-66
#>  2 infl     unemp.l1 -0.281      0.0622    -4.52  1.05e- 5
#>  3 infl     rate.l1  -0.0571     0.0295    -1.93  5.45e- 2
#>  4 infl     infl.l2  -0.542      0.0613    -8.84  4.41e-16
#>  5 infl     unemp.l2  0.234      0.0592     3.95  1.09e- 4
#>  6 infl     rate.l2   0.0580     0.0291     2.00  4.72e- 2
#>  7 infl     const     0.287      0.0955     3.01  2.96e- 3
#>  8 unemp    infl.l1  -0.00952    0.0643    -0.148 8.82e- 1
#>  9 unemp    unemp.l1  1.33       0.0670    19.8   1.09e-49
#> 10 unemp    rate.l1  -0.0210     0.0318    -0.660 5.10e- 1
#> # ℹ 11 more rows
tidy(fit, conf.int = TRUE)
#> # A tibble: 21 × 8
#>    response term     estimate std.error statistic  p.value  conf.low conf.high
#>    <chr>    <chr>       <dbl>     <dbl>     <dbl>    <dbl>     <dbl>     <dbl>
#>  1 infl     infl.l1   1.54       0.0596    25.8   2.13e-66  1.42       1.66   
#>  2 infl     unemp.l1 -0.281      0.0622    -4.52  1.05e- 5 -0.403     -0.158  
#>  3 infl     rate.l1  -0.0571     0.0295    -1.93  5.45e- 2 -0.115      0.00112
#>  4 infl     infl.l2  -0.542      0.0613    -8.84  4.41e-16 -0.663     -0.421  
#>  5 infl     unemp.l2  0.234      0.0592     3.95  1.09e- 4  0.117      0.350  
#>  6 infl     rate.l2   0.0580     0.0291     2.00  4.72e- 2  0.000715   0.115  
#>  7 infl     const     0.287      0.0955     3.01  2.96e- 3  0.0989     0.475  
#>  8 unemp    infl.l1  -0.00952    0.0643    -0.148 8.82e- 1 -0.136      0.117  
#>  9 unemp    unemp.l1  1.33       0.0670    19.8   1.09e-49  1.20       1.46   
#> 10 unemp    rate.l1  -0.0210     0.0318    -0.660 5.10e- 1 -0.0838     0.0418 
#> # ℹ 11 more rows
glance(fit)
#> # A tibble: 1 × 8
#>       K     p  nobs logLik   AIC   BIC    HQ max_root
#>   <int> <int> <int>  <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1     3     2   213  -338.  731.  822.  768.    0.946
augment(fit)
#> # A tibble: 215 × 9
#>     infl unemp  rate .fitted_infl .resid_infl .fitted_unemp .resid_unemp
#>    <dbl> <dbl> <dbl>        <dbl>       <dbl>         <dbl>        <dbl>
#>  1 1.69    2.9  1.96       NA         NA              NA         NA     
#>  2 1.45    2.7  2.19       NA         NA              NA         NA     
#>  3 1.01    2.6  2.04        1.51      -0.498           2.80      -0.198 
#>  4 0.981   3.1  1.38        0.969      0.0114          2.76       0.339 
#>  5 1.33    4.9  1.18        1.02       0.309           3.46       1.44  
#>  6 1.30    5.9  0.97        1.17       0.131           5.61       0.290 
#>  7 0.909   5.8  0.72        1.06      -0.151           6.17      -0.375 
#>  8 0.709   5.7  0.98        0.747     -0.0373          5.61       0.0890
#>  9 0.685   4.9  1.22        0.625      0.0593          5.49      -0.593 
#> 10 1.06    4.7  1.59        0.898      0.165           4.47       0.228 
#> # ℹ 205 more rows
#> # ℹ 2 more variables: .fitted_rate <dbl>, .resid_rate <dbl>
```
