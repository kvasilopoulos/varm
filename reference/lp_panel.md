# Panel local projections with fixed effects and clustered errors

For a panel \\(i,t)\\, estimates for each horizon \$\$y\_{i,t+h} =
\alpha\_{i,h} + \delta\_{t,h} + \beta_h s\_{it} +
\sum\_{l=1}^{L}(\phi\_{l,h} y\_{i,t-l} + \psi\_{l,h}s\_{i,t-l}) +
\kappa_h' x\_{it} + e\_{i,t+h}\$\$ by within (demeaning) estimation with
unit effects \\\alpha_i\\, time effects \\\delta_t\\, both, or neither.
With two-way effects on an unbalanced horizon sample the demeaning is
iterated (alternating projections), a single pass for balanced panels,
so no dummy matrix is ever built. Standard errors cluster on the unit
(or time) with the Stata-type small-sample correction
\\\frac{G}{G-1}\frac{n-1}{n-k}\\, \\k\\ counting the regressors and the
absorbed effects (identical to `lm` with dummies and hand-coded cluster
scores). The joint cross-horizon cluster covariance (scores of a cluster
summed over time, stacked over horizons) is returned for
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md).

## Usage

``` r
lp_panel(
  data,
  y,
  shock,
  unit,
  time,
  horizon = 8L,
  lags = 1L,
  controls = NULL,
  fe = c("twoway", "unit", "time", "none"),
  cluster = c("unit", "time"),
  cumulative = FALSE,
  level = c(0.68, 0.9)
)
```

## Arguments

- data:

  a data frame in long format.

- y, shock, unit, time:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  the columns of the outcome, the shock, the unit id and the time index
  (regularly spaced).

- horizon:

  maximum horizon.

- lags:

  number of lags of `y` and `shock` among the controls.

- controls:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  further (contemporaneous) control columns.

- fe:

  `"twoway"`, `"unit"`, `"time"` or `"none"` (common intercept).

- cluster:

  `"unit"` or `"time"`.

- cumulative:

  cumulate the outcome over the horizon.

- level:

  confidence levels.

## Value

an `lp_panel` object with `$estimate`, `$se`, `$bands` (`(H+1) x 1`),
`$vcov`, `$nobs`, `$nclusters`.

## Details

**Economic intuition.** Panels (countries, states, firms) give
identification and power: unit effects remove time-invariant
heterogeneity and time effects remove common shocks, so the response is
identified from cross-sectional differences in exposure to \\s\\.
Clustering by unit allows arbitrary serial correlation within the unit,
which is essential because the LP residual is an MA(\\h\\) process.
*Pitfalls:* fixed effects with lagged dependent variables induce Nickell
bias when \\T\\ is small (use long panels or drop `y` lags); time
effects absorb the average response, so the estimate is a *relative*
effect; with few clusters (\< 30) the cluster-robust bands are too
narrow; heterogeneous responses are averaged with variance weights.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
set.seed(1)
d <- expand.grid(id = 1:20, t = 1:40)
d$s <- rnorm(nrow(d))
d$y <- 0.5 * d$s + rnorm(nrow(d))
lp_panel(d, y, s, id, t, horizon = 3, lags = 1)$estimate
#>              y
#> 0  0.550833521
#> 1  0.076868592
#> 2 -0.006838165
#> 3  0.036043860
```
