# Response of one variable to one shock at a fixed horizon, across all dates

Evolution over time of \\\Theta\_{h}^{ij}(t)\\, the response of variable
\\i\\ to recursive shock \\j\\ at horizon \\h\\ computed with frozen
date-\\t\\ parameters (see
[`irf.tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md)).
This is the usual "time-varying IRF" plot: for example the impact (or
8-quarter) response of inflation to a policy shock through time. With
`normalise_to` the shock is rescaled each date so the impact response of
`variable` equals `size` (e.g. 1 percentage point on the policy rate),
removing changes in shock size.

## Usage

``` r
irf_path(
  object,
  shock,
  response,
  horizon = 0L,
  normalise_to = NULL,
  size = 1,
  level = 0.68
)
```

## Arguments

- object:

  a
  [`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
  fit.

- shock, response:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  one variable each.

- horizon:

  horizon of the response (0 = impact).

- normalise_to:

  optional
  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variable whose impact response is normalised to `size`.

- size:

  normalised impact (default 1).

- level:

  credible level of the band.

## Value

a `tvpvar_path` tibble with columns `date`, `estimate` (posterior
median), `conf.low`, `conf.high`;
[autoplot()](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
plots it.

## See also

Other TVP-VAR:
[`tidy.tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md),
[`tvp_dynamics`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md),
[`tvp_prior()`](https://kvasilopoulos.github.io/varm/reference/tvp_prior.md),
[`tvp_volatility()`](https://kvasilopoulos.github.io/varm/reference/tvp_volatility.md),
[`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)

## Examples

``` r
# \donttest{
fit <- tvpvar(yraw, p = 2, ndraw = 200, nburn = 200)
irf_path(fit, shock = rate, response = infl, horizon = 4)
#> # A tibble: 173 × 4
#>     date estimate conf.low conf.high
#>  * <int>    <dbl>    <dbl>     <dbl>
#>  1    43  -0.0126  -0.0375   0.0101 
#>  2    44  -0.0126  -0.0354   0.0105 
#>  3    45  -0.0124  -0.0330   0.0102 
#>  4    46  -0.0126  -0.0355   0.00945
#>  5    47  -0.0136  -0.0353   0.00917
#>  6    48  -0.0126  -0.0353   0.00987
#>  7    49  -0.0123  -0.0348   0.0106 
#>  8    50  -0.0124  -0.0348   0.0106 
#>  9    51  -0.0127  -0.0380   0.0105 
#> 10    52  -0.0123  -0.0423   0.0116 
#> # ℹ 163 more rows
# }
```
