# State-dependent local projections (Auerbach-Gorodnichenko)

Local projections whose coefficients depend on the state of the economy
through a transition function of a switching variable \\z\\. For each
horizon \\h\\, \$\$y\_{t+h} = F(z\_{t-1})\\\[\alpha^L_h + \beta^L_h
s_t + \gamma^{L\prime}\_h w_t\] + (1-F(z\_{t-1}))\\\[\alpha^H_h +
\beta^H_h s_t + \gamma^{H\prime}\_h w_t\] + e\_{t+h},\$\$ with the
smooth logistic transition \\F(z) = 1/(1+\exp\\\gamma(\tilde z - c)\\)\\
on the standardised variable \\\tilde z=(z-\bar z)/\mathrm{sd}(z)\\
(Auerbach & Gorodnichenko 2012, 2013; \\F\to1\\ in the "low" state, e.g.
recession when \\z\\ is growth) or the binary version \\F = 1\\\tilde
z\_{t-1} \< c\\\\. The state enters lagged by one period to avoid
contemporaneous feedback of the shock on the state. All \\2\times\\
regressors are estimated in one regression per horizon (C++ engine),
with Newey-West errors, and the cross-horizon covariance is used to test
\\H_0:\beta^L_h=\beta^H_h\\ pointwise (HAC t-test) and jointly across
horizons (Wald).

## Usage

``` r
lp_state(
  Y,
  shock,
  state,
  horizon = 12L,
  lags = 4L,
  contemp = NULL,
  controls = NULL,
  responses = NULL,
  cumulative = NULL,
  transition = c("logistic", "binary"),
  gamma = 1.5,
  cutoff = 0,
  low_share = NULL,
  hp_lambda = NULL,
  standardise = TRUE,
  timing = c("contemporaneous", "lagged"),
  lags_lin = lags,
  shock_type = 0L,
  state_const = TRUE,
  trend = FALSE,
  se = c("nw", "hc"),
  nw_lag = NULL,
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

- state:

  switching variable \\z\\: a variable of `Y` (tidyselect, one variable)
  or a numeric vector of length T.

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

- responses:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  responses to use (default all variables of `Y`).

- cumulative:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variables whose response is cumulated (\\\sum\_{j=0}^h y\_{t+j}\\),
  for variables measured in differences.

- transition:

  `"logistic"` (smooth) or `"binary"`.

- gamma:

  logistic smoothness (larger is closer to a regime switch).

- cutoff:

  cutoff \\c\\ on the standardised scale; ignored if `low_share` is
  given.

- low_share:

  optional share in (0, 1) of observations to place below the cutoff
  (sets \\c\\ to that quantile of \\\tilde z\\).

- hp_lambda:

  if not `NULL`, use the Hodrick-Prescott cycle of `state` (smoothing
  parameter, e.g. 1600 for quarterly data) as switching variable.

- standardise:

  standardise \\z\\ (after filtering) to mean 0, sd 1.

- timing:

  `"contemporaneous"` or `"lagged"` (the `lpirfs` convention), see
  Details.

- lags_lin:

  lag length of the linear VAR used for the impact matrix
  (`timing = "lagged"`).

- shock_type:

  `"lagged"` only: `0` standard-deviation, `1` unit shock.

- state_const:

  interact the constant with the state (state-specific intercepts);
  `lpirfs` uses a common constant, so set `FALSE` when comparing.

- trend:

  `"lagged"` only: include a linear time trend (`lpirfs` trend = 1).

- se:

  `"nw"` Newey-West/Bartlett or `"hc"` Eicker-Huber-White (HC1).

- nw_lag:

  Newey-West truncation lag; default `h + 1`.

- level:

  confidence levels for the pointwise bands.

- dof_adjust:

  apply an `n/(n-m)` small-sample factor to the Newey-West variance.

## Value

an `lp_state` object: arrays `$estimate`, `$se` (`(H+1) x K x 3`, last
dimension `high`, `low`, `diff = low - high`), `$bands`,
`$tstat`/`$pvalue` for the difference, `$wald` (joint test of equality
over horizons, per response), the transition `$F` and `$cutoff`, and
`$raw` (coefficient array and cross-horizon covariance from the engine)
with `$contrasts`.
[tidy()](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
returns one row per horizon, response, state and band level.

## Details

**Timing conventions.** `timing = "contemporaneous"` (default) is the
shock-regression form above (observed or external shock \\s_t\\, as in
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md)).
`timing = "lagged"` replicates the convention of
[`lpirfs::lp_nl`](https://rdrr.io/pkg/lpirfs/man/lp_nl.html): the
regressors are lags of \\y\\ interacted with the state, with no
contemporaneous terms; the horizon-\\h\\ response is the coefficient on
\\y\_{t-1}\\ in the equation for \\y\_{t+h-1}\\ multiplied by column
`shock` of the Cholesky impact matrix of a linear VAR
(standard-deviation shock when `shock_type = 0`, unit own-impact when
`1`), and the impact response (h = 0) is that column of the impact
matrix in both states. Differences to `lpirfs`: bands here are normal
quantiles at `level` (not a user-supplied multiplier); `lpirfs` uses
Newey-West lag \\h\\ for its horizon \\h\\ (regression of
\\y\_{t+h-1}\\), which equals the default \\L=(h-1)+1\\ here; `lpirfs`
has no state-specific constants (`state_const = FALSE`).

**Economic intuition.** Fiscal and monetary multipliers are often larger
in recessions; this estimates state-specific responses without a
nonlinear VAR, so the state can change over the horizon without being
imposed. The cost is variance: each state uses only the observations
weighted by \\F\\ or \\1-F\\, so effective samples are small and bands
wide, and with persistent states the difference test has low power at
long horizons. *Pitfalls:* the calibrated `gamma` and cutoff are not
estimated (AG use \\\gamma=1.5\\ and a cutoff so that the low state has
about 20% weight, see `low_share`); a state that is endogenous to the
shock beyond one lag contaminates the comparison; with the binary
version one state may have very few observations.

## References

Auerbach, A. and Gorodnichenko, Y. (2012) Measuring the output responses
to fiscal policy. *AEJ: Economic Policy*; (2013) Fiscal multipliers in
recession and expansion.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
fit <- lp_state(
  yraw,
  shock = rate,
  state = infl,
  horizon = 8,
  lags = 2,
  contemp = c(infl, unemp),
  gamma = 3,
  low_share = 0.3
)
fit$estimate[, 3, ]
#>         high       low         diff
#> 0 1.00000000 1.0000000 2.831069e-13
#> 1 0.58913461 1.6408691 1.051734e+00
#> 2 0.32338423 1.8500586 1.526674e+00
#> 3 0.37654588 1.6185735 1.242028e+00
#> 4 0.45306705 1.4506486 9.975816e-01
#> 5 0.62270301 1.1966626 5.739596e-01
#> 6 0.42389376 0.8809570 4.570632e-01
#> 7 0.12435717 0.6974001 5.730430e-01
#> 8 0.09748667 0.4688301 3.713435e-01
```
