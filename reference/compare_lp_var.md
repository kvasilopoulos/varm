# Compare local-projection and VAR impulse responses

Plagborg-Moller & Wolf (2021): with the same lag length and the same
recursive identification, population LP and VAR responses coincide at
all horizons; in finite samples they differ only by estimation error.
This helper estimates the LP
\\y\_{t+h}=\alpha_h+\beta_hy\_{j,t}+\gamma_h'w_t+e\\ with controls the
contemporaneous variables ordered before the shock \\j\\ and \\p\\ lags
of all variables, and the VAR(\\p\\) Cholesky response (ordering =
column order) normalised to a unit impact on the shock variable,
\\\Theta_h e_j/\Theta_0^{jj}\\, and returns both with their difference.

## Usage

``` r
compare_lp_var(
  Y,
  shock,
  horizon = 20L,
  lags = 4L,
  se = c("nw", "hc"),
  nw_lag = NULL
)
```

## Arguments

- Y:

  data frame / tibble / matrix (T x K) of responses (and the source of
  lag controls). Non-numeric columns (e.g. a date) are dropped.

- shock:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  the shock variable of `Y`.

- horizon:

  maximum horizon H.

- lags:

  number of lags of `Y` in the controls (and of an external shock).

- se:

  `"nw"` Newey-West/Bartlett or `"hc"` Eicker-Huber-White (HC1).

- nw_lag:

  Newey-West truncation lag; default `h + 1`.

## Value

a `lp_var_cmp` tibble with `horizon`, `response`, `lp`, `lp_se`, `var`,
`diff = lp - var`, `z = diff / lp_se`.

## Details

**Economic intuition.** Choosing between LP and VAR is a bias-variance
trade-off: the VAR extrapolates the first \\p\\ lags to all horizons
(low variance, bias if the lag order is too short for long-run
dynamics), LP is unbiased for any \\h\\ but has variance growing with
\\h\\. A divergence at long horizons that is large relative to `lp_se`
signals VAR misspecification (or LP noise); agreement supports either.
*Pitfall:* the comparison is only meaningful with identical lags and
ordering, which this function enforces. VAR responses are point
estimates (no bands).

## References

Plagborg-Moller, M. and Wolf, C. (2021) Local projections and VARs
estimate the same impulse responses. *Econometrica*.

## See also

Other local projections:
[`autoplot.lpm()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md),
[`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md),
[`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md),
[`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md),
[`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md),
[`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md),
[`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md),
[`tidy.lpm()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)

## Examples

``` r
compare_lp_var(yraw, shock = unemp, horizon = 12, lags = 2)
#> <lp_var_cmp> shock 'unemp', 2 lags; mean |LP - VAR| / se_LP by horizon block:
#>         infl  rate unemp
#> short  0.252 0.276 0.377
#> medium 0.498 0.152 0.344
#> long   0.189 0.519 0.500
```
