# MCMC convergence diagnostics

Effective sample size (Geyer's initial monotone positive sequence,
multi-chain form of Vehtari et al. 2021), split R-hat, and the Geweke
(1992) mean-equality z-score (first 10% vs last 50% of each chain,
variances inflated by the autocorrelation time) for every parameter.

## Usage

``` r
mcmc_diagnostics(draws)

# S3 method for class 'mcmc_diagnostics'
tidy(x, ...)

# S3 method for class 'mcmc_diagnostics'
as_tibble(x, ...)

# S3 method for class 'mcmc_diagnostics'
glance(x, ...)
```

## Arguments

- draws:

  one of: a numeric matrix `iterations x parameters` (one chain); a 3-d
  array `iterations x chains x parameters`; a list of such matrices (one
  per chain); a `bvarm` (hyperparameter draws if hierarchical, otherwise
  the coefficient draws `B`); or a `varm_draws` object (coefficients).

- x:

  an `mcmc_diagnostics` object.

- ...:

  unused.

## Value

a tibble of class `mcmc_diagnostics` with columns `term`, `estimate`
(posterior mean), `std.error` (posterior sd), `ess` (bulk effective
sample size over all chains), `rhat` (split R-hat; every chain is cut in
two halves, so a single chain works) and `geweke_z` (mean over chains,
`|z| > 2` suggests non-stationarity). The draws are kept for
[autoplot()](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.md);
[`glance()`](https://generics.r-lib.org/reference/glance.html)
summarises the table.

## See also

Other Bayesian VAR:
[`autoplot.mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.md),
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md),
[`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.md),
[`prior_triangular()`](https://kvasilopoulos.github.io/varm/reference/prior_triangular.md),
[`tidy.bvarm()`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md),
[`volatility()`](https://kvasilopoulos.github.io/varm/reference/volatility.md)

## Examples

``` r
set.seed(1)
x <- cbind(a = arima.sim(list(ar = 0.8), 2000), b = rnorm(2000))
diag <- mcmc_diagnostics(x)
diag
#>  term estimate std.error  ess rhat geweke_z
#>     a  -0.0902      1.72  219    1    0.603
#>     b   0.0174      1.03 1865    1   -0.569
glance(diag)
#> # A tibble: 1 × 6
#>   n_param n_chain n_iter min_ess max_rhat max_abs_geweke
#>     <int>   <int>  <int>   <dbl>    <dbl>          <dbl>
#> 1       2       1   2000    219.    1.000          0.603
```
