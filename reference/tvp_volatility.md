# Volatility paths of a TVP-VAR

Posterior bands for the time-varying standard deviations:
`type = "reduced"` gives \\\sqrt{\Sigma\_{t,ii}}\\ (reduced-form
residual s.d. of each variable, what a forecaster faces); `"structural"`
gives \\\exp(h\_{it}/2)\\, the s.d. of the orthogonalised (recursive)
shocks. Comparing the two separates changes in shock size from changes
in contemporaneous propagation \\A_t\\.

## Usage

``` r
tvp_volatility(object, type = c("reduced", "structural"), level = 0.68)
```

## Arguments

- object:

  a
  [`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
  fit.

- type:

  `"reduced"` or `"structural"`.

- level:

  credible level.

## Value

a tibble with `date`, `variable`, `estimate` (posterior median),
`conf.low`, `conf.high`.

## See also

Other TVP-VAR:
[`irf_path()`](https://kvasilopoulos.github.io/varm/reference/irf_path.md),
[`tidy.tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md),
[`tvp_dynamics`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md),
[`tvp_prior()`](https://kvasilopoulos.github.io/varm/reference/tvp_prior.md),
[`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
