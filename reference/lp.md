# Local projections

Jorda (2005) local projections: for each horizon \\h=0,\dots,H\\ and
each response \\y_k\\, \$\$y\_{k,t+h} = \alpha_h + \beta_h s_t +
\gamma_h' w_t + e\_{t+h},\$\$ where \\s_t\\ is the shock, \\w_t\\ the
controls (lags of `Y`, optional contemporaneous variables, extra
controls). \\\beta_h\\ is the impulse response at horizon \\h\\. The
regressions, 2SLS and HAC standard errors run in C++.

## Usage

``` r
lp(
  Y,
  shock,
  horizon = 20L,
  lags = 4L,
  contemp = NULL,
  controls = NULL,
  iv = NULL,
  responses = NULL,
  cumulative = NULL,
  se = c("nw", "hc"),
  nw_lag = NULL,
  lag_augment = FALSE,
  level = c(0.68, 0.9),
  dof_adjust = FALSE
)
```

## Arguments

- Y:

  data frame / tibble / matrix (T x K) of responses (and the source of
  lag controls). Non-numeric columns (e.g. a date) are dropped.

- shock:

  the shock \\s_t\\: a variable of `Y` (tidyselect; observed shock,
  unit-impact normalisation) or a numeric vector of length T (external
  shock series, e.g. narrative).

- horizon:

  maximum horizon H.

- lags:

  number of lags of `Y` in the controls (and of an external shock).

- contemp:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variables of `Y` entered contemporaneously as controls. For a
  recursive (Cholesky) identification with `shock` a variable of `Y`,
  list the variables ordered *before* the shock.

- controls:

  optional extra T x q matrix of controls.

- iv:

  optional T x q matrix of instruments for `shock` (LP-IV; Stock &
  Watson 2018).

- responses:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  responses to use (default all variables of `Y`).

- cumulative:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variables whose response is cumulated (\\\sum\_{j=0}^h y\_{t+j}\\),
  for variables measured in differences.

- se:

  `"nw"` Newey-West/Bartlett or `"hc"` Eicker-Huber-White (HC1).

- nw_lag:

  Newey-West truncation lag; default `h + 1`.

- lag_augment:

  Montiel Olea & Plagborg-Moller (2021) lag augmentation: add one extra
  lag that is not interpreted and use `"hc"` standard errors, valid
  without stationarity assumptions and robust at long horizons.

- level:

  confidence levels for the pointwise bands.

- dof_adjust:

  apply an `n/(n-m)` small-sample factor to the Newey-West variance.

## Value

an `lpm` object: `$estimate` (`(H+1) x K`), `$se`, `$bands`, `$nobs`.
Use [tidy()](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
for a tibble and
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
for a `ggplot`.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
fit <- lp(yraw, shock = rate, horizon = 12, lags = 2, contemp = c(infl, unemp))
fit
#> <lpm> local projection of infl, unemp, rate on 'rate', horizons 0..12 (2 lags, NW)
#>       infl   unemp    rate
#> 0   0.0000  0.0000  1.0000
#> 1  -0.0589 -0.0295  0.7650
#> 4  -0.0217  0.0284  0.5960
#> 8  -0.0745  0.1937  0.1735
#> 12 -0.4550  0.1486 -0.1578
tidy(fit)
#> # A tibble: 78 × 7
#>    horizon response  estimate std.error level  conf.low conf.high
#>      <int> <fct>        <dbl>     <dbl> <dbl>     <dbl>     <dbl>
#>  1       0 infl      1.51e-14  9.62e-14   0.9 -1.43e-13  1.73e-13
#>  2       1 infl     -5.89e- 2  4.20e- 2   0.9 -1.28e- 1  1.01e- 2
#>  3       2 infl     -6.16e- 2  7.23e- 2   0.9 -1.81e- 1  5.74e- 2
#>  4       3 infl     -3.15e- 2  1.05e- 1   0.9 -2.04e- 1  1.41e- 1
#>  5       4 infl     -2.17e- 2  1.52e- 1   0.9 -2.72e- 1  2.28e- 1
#>  6       5 infl     -7.74e- 3  1.71e- 1   0.9 -2.88e- 1  2.73e- 1
#>  7       6 infl     -1.53e- 2  1.79e- 1   0.9 -3.10e- 1  2.79e- 1
#>  8       7 infl     -4.27e- 2  1.91e- 1   0.9 -3.57e- 1  2.71e- 1
#>  9       8 infl     -7.45e- 2  1.97e- 1   0.9 -3.99e- 1  2.50e- 1
#> 10       9 infl     -1.68e- 1  1.81e- 1   0.9 -4.65e- 1  1.30e- 1
#> # ℹ 68 more rows
```
