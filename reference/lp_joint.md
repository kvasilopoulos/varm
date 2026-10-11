# Local projections with joint cross-horizon covariance

Same regressions as
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md) but
additionally returns, for every response, the covariance matrix of the
whole impulse-response vector \\(\hat\beta_0,\dots,\hat\beta_H)'\\,
\$\$\widehat{\mathrm{Cov}}(\hat\beta_h,\hat\beta\_{h'}) =
\sum\_{l=-L}^{L} w_l \sum_t q\_{t,h} e\_{t+h}\\ q\_{t-l,h'}
e\_{t-l+h'},\quad q\_{t,h} = x\_{t,h}'(X_h'X_h)^{-1}e_1,\$\$ with
Bartlett weights \\w_l = 1-\|l\|/(L+1)\\ and \\L = \max(h,h')+1\\ (or a
common `nw_lag`). The diagonal reproduces the standard errors of
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md). The
engine is `cpp_lp_cov` (C++).

## Usage

``` r
lp_joint(
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

an `lp_joint` object (also an `lpm`) with the extra element `$vcov`: a
named list (one per response) of \\(H+1)\times(H+1)\\ covariance
matrices.

## Details

**Economic intuition.** Local-projection estimates at neighbouring
horizons are built from overlapping data and highly correlated, so
reading a pointwise band as a statement about the whole response path
overstates what the data say. The joint covariance is the input for
simultaneous inference
([`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)),
joint hypothesis tests (e.g. "the response is zero at all horizons"),
and smoothing. *Pitfall:* with the default horizon-specific bandwidth
the full matrix is not guaranteed positive semi-definite; use a common
`nw_lag` if you need to invert it.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
fit <- lp_joint(yraw, shock = rate, horizon = 8, lags = 2, contemp = c(infl, unemp))
sqrt(diag(fit$vcov[[1]])) - fit$se[, 1] # zero
#>             0             1             2             3             4 
#> -2.524355e-29  2.775558e-17  1.387779e-17  4.163336e-17  5.551115e-17 
#>             5             6             7             8 
#>  2.775558e-17  2.775558e-17  0.000000e+00  0.000000e+00 
```
