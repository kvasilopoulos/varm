# Direct multi-step forecasts from local projections

For each horizon \\h=1,\dots,H\\ regresses \\y\_{t+h}\\ on
\\x_t=(1,y_t',\dots,y\_{t-p+1}')'\\ and forecasts \\\hat
y\_{T+h}=x_T'\hat B_h\\. The forecast standard error is
\\\hat\sigma_h\sqrt{1+x_T'(X_h'X_h)^{-1}x_T}\\ (the direct-regression
residual has an MA(h-1) structure but its unconditional variance is what
matters for the forecast error).

## Usage

``` r
lp_forecast(Y, horizon = 8L, lags = 4L, level = c(0.68, 0.9))
```

## Arguments

- Y:

  data frame / tibble / matrix (T x K); non-numeric columns are dropped.

- horizon:

  maximum forecast horizon.

- lags:

  number of lags \\p\\.

- level:

  confidence levels for the bands.

## Value

an `lp_forecast` object: `$mean`, `$se` (`H x K`), `$bands`
(`H x K x 2*length(level)`).
[tidy()](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
returns a tibble.

## Details

**Economic intuition.** Iterated (VAR) forecasts are efficient if the
VAR is correctly specified but compound any misspecification with the
horizon; direct forecasts estimate each horizon's conditional mean
separately, which is robust to misspecification but noisier (Marcellino,
Stock & Watson 2006). Prefer direct forecasts at long horizons with long
samples and a possibly nonlinear/mis-specified dynamics, iterated
forecasts when the sample is short. At \\h=1\\ the two coincide exactly.
*Pitfall:* forecasts at different horizons are not constrained to be
consistent with each other and can be jagged.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
lp_forecast(yraw, horizon = 4, lags = 2)
#> <lp_forecast> direct forecasts, horizons 1..4
#>     infl  unemp   rate
#> 1 2.7209 4.9214 4.6197
#> 2 2.6601 5.0377 4.4461
#> 3 2.6606 5.2093 4.2187
#> 4 2.7711 5.3189 4.2035
```
