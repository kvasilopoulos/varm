# Smooth local projections (Barnichon & Brownlees)

Imposes that the impulse response is a smooth function of the horizon,
\\\beta(h) = \sum\_{j=1}^J \theta_j B_j(h)\\ with a B-spline basis
\\B_j\\ (cubic by default, equally spaced knots over \\0..H\\), and
estimates all horizons jointly by penalised least squares:
\$\$\hat\theta = \arg\min\_\theta \sum\_{h=0}^H\sum_t \Big(\tilde
y\_{t+h}-\sum_j\theta_jB_j(h)\tilde s\_{t,h}\Big)^2 +
\lambda\\\theta'D_r'D_r\theta,\$\$ where tilde denotes residuals from
the horizon-specific regression on the controls (Frisch-Waugh), \\D_r\\
is the \\r\\-th order difference matrix (P-spline penalty, Eilers & Marx
1996) which shrinks \\\beta(h)\\ towards a polynomial of degree \\r-1\\
as \\\lambda\to\infty\\. \\\lambda\\ is chosen by K-fold
cross-validation over contiguous blocks of time (preserving serial
dependence). Standard errors use the HAC sandwich \\(A+\lambda
P)^{-1}\hat\Omega(A+\lambda P)^{-1}\\, \\\hat\Omega\\ Newey-West on the
scores summed over horizons at each \\t\\; they treat the smoothing bias
as negligible.

## Usage

``` r
lp_smooth(
  Y,
  shock,
  response = 1L,
  horizon = 20L,
  lags = 4L,
  contemp = NULL,
  controls = NULL,
  cumulative = FALSE,
  nbasis = NULL,
  degree = 3L,
  penalty_order = 2L,
  lambda = NULL,
  nfolds = 5L,
  nw_lag = NULL,
  level = c(0.68, 0.9)
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

- response:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  the single response variable of `Y`.

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

- cumulative:

  cumulate the response (for variables in differences).

- nbasis:

  number of B-spline basis functions \\J\\ (default `min(H + 1, 8)`).

- degree:

  spline degree.

- penalty_order:

  order \\r\\ of the difference penalty.

- lambda:

  a fixed (absolute) penalty; `NULL` = choose by cross-validation over a
  grid scaled by \\\mathrm{tr}(A)/J\\.

- nfolds:

  number of CV blocks.

- nw_lag:

  Newey-West truncation for the score (default `H + 1`).

- level:

  confidence levels for the pointwise bands.

## Value

an `lp_smooth` object: `$estimate`/`$se` (`(H+1) x 1`), `$bands`,
`$vcov` (compatible with
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)),
`$lp`/`$lp_se` (unsmoothed LP), `$lambda`, `$cv` (grid and CV errors),
`$edf`.

## Details

**Economic intuition.** Raw LP estimates are noisy and jagged because
every horizon is estimated on its own; theory says responses are smooth
(hump-shaped, monotone decay). Smoothing borrows strength across
neighbouring horizons, trading a small bias (a hump is flattened if
\\\lambda\\ is too large) for a large variance reduction, which is what
the LP-vs-VAR bias-variance trade-off is about: smooth LP sits between
unrestricted LP and a parametric VAR. *Pitfalls:* sharp features (e.g.
an impact jump at \\h=0\\, a sudden reversal) are smoothed away; use
enough basis functions; the bands are conditional on \\\lambda\\,
ignoring the CV step; one response at a time.

## References

Barnichon, R. and Brownlees, C. (2019) Impulse response estimation by
smooth local projections. *Review of Economics and Statistics*.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
fit <- lp_smooth(
  yraw,
  shock = rate,
  response = infl,
  horizon = 12,
  lags = 2,
  contemp = c(infl, unemp)
)
round(cbind(smooth = fit$estimate[, 1], raw = fit$lp), 3)
#>    smooth    raw
#> 0   0.011  0.000
#> 1   0.003 -0.059
#> 2  -0.007 -0.062
#> 3  -0.019 -0.032
#> 4  -0.034 -0.022
#> 5  -0.054 -0.008
#> 6  -0.080 -0.015
#> 7  -0.114 -0.043
#> 8  -0.155 -0.075
#> 9  -0.203 -0.168
#> 10 -0.256 -0.276
#> 11 -0.313 -0.378
#> 12 -0.372 -0.455
```
