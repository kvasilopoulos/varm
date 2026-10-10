# Posterior volatility paths

Posterior mean of the stochastic volatilities \\\exp(h\_{it}/2)\\
(standard deviation of the structural shock of equation i) from a
`bvarm` fitted with `sv = TRUE`.

## Usage

``` r
volatility(object)
```

## Arguments

- object:

  a `bvarm` with stochastic volatility.

## Value

`T_eff x K` matrix.

## See also

Other Bayesian VAR:
[`autoplot.mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.md),
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md),
[`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md),
[`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.md),
[`prior_triangular()`](https://kvasilopoulos.github.io/varm/reference/prior_triangular.md),
[`tidy.bvarm()`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md)
