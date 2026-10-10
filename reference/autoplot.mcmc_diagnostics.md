# Plot MCMC diagnostics

`type = "summary"` shows the effective sample size and R-hat of every
parameter, `"trace"` the chains and `"density"` the marginal posterior
densities of the selected parameters.

## Usage

``` r
# S3 method for class 'mcmc_diagnostics'
autoplot(
  object,
  type = c("summary", "trace", "density"),
  params = NULL,
  max_params = 12L,
  ...
)

# S3 method for class 'mcmc_diagnostics'
plot(x, ...)
```

## Arguments

- object:

  an
  [`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md)
  result.

- type:

  `"summary"`, `"trace"` or `"density"`.

- params:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  parameters for the trace and density plots (default: the first
  `max_params`).

- max_params:

  number of parameters shown when `params` is `NULL`.

- ...:

  unused.

- x:

  an `mcmc_diagnostics` object.

## Value

a `ggplot` object.

## See also

Other Bayesian VAR:
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md),
[`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md),
[`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.md),
[`prior_triangular()`](https://kvasilopoulos.github.io/varm/reference/prior_triangular.md),
[`tidy.bvarm()`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md),
[`volatility()`](https://kvasilopoulos.github.io/varm/reference/volatility.md)

## Examples

``` r
set.seed(1)
d <- mcmc_diagnostics(cbind(a = rnorm(500), b = rnorm(500)))
autoplot(d)

autoplot(d, "trace")
```
