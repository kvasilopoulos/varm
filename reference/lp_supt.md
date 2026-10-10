# Sup-t simultaneous confidence bands for impulse responses

Montiel Olea & Plagborg-Moller (2019) sup-t bands. With \\\hat\sigma_h\\
the standard error of \\\hat\beta_h\\ and \\\hat R\\ the correlation
matrix implied by the joint covariance, the band is \$\$\hat\beta_h \pm
\hat c\_{1-\alpha}\hat\sigma_h,\qquad \hat c\_{1-\alpha} =
(1-\alpha)\text{-quantile of } \max_h \|z_h\|,\\ z\sim N(0,\hat R),\$\$
computed by simulation. The band covers the entire response path with
probability \\1-\alpha\\ (asymptotically), is narrower than Bonferroni
and is never narrower than the pointwise band.

## Usage

``` r
lp_supt(object, level = 0.9, nsim = 10000L, seed = NULL)
```

## Arguments

- object:

  an object with a `$vcov` element:
  [`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
  [`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md)
  or
  [`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md).

- level:

  simultaneous coverage level.

- nsim:

  number of simulation draws.

- seed:

  optional RNG seed.

## Value

a tibble with columns `response`, `horizon`, `estimate`, `std.error`,
`pw_low`, `pw_high` (pointwise at `level`), `conf.low`, `conf.high`
(sup-t); the critical values are in `attr(, "crit")`.
[tidy()](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
stacks the two band types and
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
compares them.

## Details

**Economic intuition.** A pointwise 90% band is a statement about one
horizon at a time; across 20 horizons the chance that the true path
leaves it somewhere is far above 10%. Use sup-t bands when you want to
claim that a whole shape (a hump, a sign, "no response at any horizon")
is supported. Because LP estimates are strongly positively correlated
across horizons the sup-t critical value is usually only modestly above
1.64 (typically 2.0-2.5), not the 3+ of Bonferroni. *Pitfall:* it needs
a good joint covariance (heavy persistence makes HAC covariances
unreliable at long horizons) and horizons with zero variance (normalised
impact) are ignored.

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
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
fit <- lp_joint(yraw, shock = rate, horizon = 8, lags = 2, contemp = c(infl, unemp))
b <- lp_supt(fit, level = 0.9, seed = 1)
attr(b, "crit")
#>     infl    unemp     rate 
#> 2.251103 2.323445 2.267141 
```
