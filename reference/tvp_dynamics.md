# Time-varying impulse responses, variance decompositions and forecasts of a TVP-VAR

Impulse responses are computed *at a date* by freezing the parameters at
their date-`t` values (the standard "local" IRF of Cogley-Sargent /
Primiceri): \$\$\Theta_h(t)=\Phi_h(B_t)\\A_t^{-1}\Lambda_t^{1/2},\qquad
\Phi_h=\sum\_{l=1}^{\min(h,p)}A_l(t)\Phi\_{h-l},\$\$ i.e. the response
of \\y\_{t+h}\\ to a one-standard-deviation recursive shock if the
economy were to stay at the date-\\t\\ parameters. Note that the shock
*size* \\\Lambda_t^{1/2}\\ varies with the date, so changes in the IRF
mix changes in propagation and in volatility; use
`normalise = list(shock =, variable =, size =)` to hold the impact
response fixed (e.g. a 25bp rate shock) and isolate changes in
propagation. The FEVD uses the same frozen-parameter responses. The
identification is recursive in the order of the data.

## Usage

``` r
# S3 method for class 'tvpvar'
irf(
  object,
  date = NULL,
  horizon = 20L,
  level = c(0.68, 0.9),
  normalise = NULL,
  cumulative = NULL,
  ...,
  .fevd = FALSE
)

# S3 method for class 'tvpvar'
fevd(object, date = NULL, horizon = 20L, level = c(0.68, 0.9), ...)

# S3 method for class 'tvpvar'
forecast(object, h = 8L, level = c(0.68, 0.9), drift = TRUE, ...)
```

## Arguments

- object:

  a
  [`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
  fit.

- date:

  one or several dates (matched to the nearest element of
  `object$dates`: the `ts` time stamp if `Y` was a `ts`, else the row
  index of `Y`), or `"last"`/`"first"`. Default: the last date.

- horizon:

  maximum horizon (impact is horizon 0).

- level:

  credible/confidence levels (e.g. `c(0.68, 0.90)`); bands are
  equal-tailed percentiles.

- normalise:

  optional
  [`normalise_to()`](https://kvasilopoulos.github.io/varm/reference/normalise_to.md)
  specification (the older `list(shock =, variable =, size =)` still
  works): rescale each shock so that the impact response of `variable`
  equals `size` (e.g. a 100bp policy-rate shock).

- cumulative:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variables whose responses are cumulated over the horizon (use for
  variables entered in differences), e.g. `cumulative = c(infl, unemp)`.

- ...:

  unused.

- .fevd:

  internal; compute the variance decomposition instead of the IRF.

- h:

  forecast horizon (steps).

- drift:

  simulate future parameter drift (`TRUE`) or freeze parameters at the
  last date.

## Value

for `irf`/`fevd`: a `varm_irf`/`varm_fevd` (single date, plottable) or a
named list of them of class `tvpvar_irf` (several dates). `forecast`
returns a `varm_forecast`.

## Details

Forecasts simulate, for every stored posterior draw, future paths of
\\\beta,\alpha,\log\lambda\\ as random walks with the draw's own
\\Q,S,W\\ (`drift = TRUE`), so predictive densities include the
possibility that parameters keep moving; `drift = FALSE` freezes them at
the end-of-sample values.

**Economic intuition.** Comparing `irf(fit, date = c(1975, 2005))` asks
whether, say, the output response to a monetary tightening was larger in
the 1970s than in the Great Moderation. Bands reflect posterior
uncertainty about the date-specific parameters only; they do not reflect
uncertainty about how long the parameters stay put.

## See also

Other TVP-VAR:
[`irf_path()`](https://kvasilopoulos.github.io/varm/reference/irf_path.md),
[`tidy.tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md),
[`tvp_prior()`](https://kvasilopoulos.github.io/varm/reference/tvp_prior.md),
[`tvp_volatility()`](https://kvasilopoulos.github.io/varm/reference/tvp_volatility.md),
[`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
