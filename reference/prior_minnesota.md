# Minnesota / sum-of-coefficients / single-unit-root prior

Prior of Giannone, Lenza & Primiceri (2015) for a conjugate
Normal-inverse-Wishart BVAR: \$\$\beta\|\Sigma \sim
N(b,\Sigma\otimes\Omega),\quad \Sigma\sim IW(\Psi, M+2).\$\$ The prior
variance of the coefficient on lag \\l\\ of variable \\j\\ in equation
\\i\\ is \\\lambda^2/(l^\alpha\psi_j)\\\psi_i\\-scaled through
\\\Sigma\\; the constant is (nearly) flat. Dummy-observation priors
shrink towards unit roots (`soc`, Doan-Litterman-Sims) and towards
cointegration/no-trend (`sur`, Sims).

## Usage

``` r
prior_minnesota(
  lambda = c(0.2, 0.4, 1e-04, 5),
  alpha = c(2, 0.25, 1, 3),
  soc = c(1, 1, 1e-04, 50),
  sur = c(1, 1, 1e-04, 50),
  psi = NULL,
  delta = 1,
  hyper = c("lambda", "soc", "sur"),
  var = 1e+07
)
```

## Arguments

- lambda:

  overall tightness (`mode` 0.2).

- alpha:

  lag-decay exponent (`mode` 2; fixed unless in `hyper`).

- soc:

  sum-of-coefficients tightness \\\mu\\ (smaller = tighter), or `NULL`
  to drop.

- sur:

  single-unit-root / dummy-initial-observation tightness \\\theta\\, or
  `NULL` to drop.

- psi:

  `NULL` (use AR(p) residual variances) or a numeric vector of length K
  (the scale of the inverse-Wishart prior, variance units).

- delta:

  prior mean of the own first lag: `1` (random walk, levels) or `0`
  (stationary, growth rates); scalar or length K.

- hyper:

  character vector among `"lambda"`, `"alpha"`, `"soc"`, `"sur"`,
  `"psi"`; or `"full"`. Default `c("lambda", "soc", "sur")`;
  `character(0)` fixes everything.

- var:

  prior variance of deterministic-term coefficients.

## Value

a `prior_minnesota` object for
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md).

## Details

Each hyperparameter is a `c(mode, sd, min, max)` vector (gamma
hyperprior with that mode and standard deviation; hard bounds). `hyper`
names the hyperparameters estimated hierarchically; all others stay
fixed at their `mode`.

## See also

Other Bayesian VAR:
[`autoplot.mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.md),
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md),
[`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md),
[`prior_triangular()`](https://kvasilopoulos.github.io/varm/reference/prior_triangular.md),
[`tidy.bvarm()`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md),
[`volatility()`](https://kvasilopoulos.github.io/varm/reference/volatility.md)
