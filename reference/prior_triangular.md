# Triangular-sampler priors (large BVAR, shrinkage, stochastic volatility)

Prior for the equation-by-equation sampler of
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
(Carriero, Clark & Marcellino 2019 style): structural equation \\i\\
regresses \\y_i\\ on lags, deterministics and the contemporaneous
\\y_1,\dots,y\_{i-1}\\. Variables are therefore *ordered*: the
reduced-form fit is order-dependent in the prior (not in the
likelihood).

## Usage

``` r
prior_triangular(
  type = c("minnesota", "ssvs", "horseshoe"),
  lambda1 = 0.2,
  lambda2 = 0.5,
  alpha = 2,
  lambda3 = 100,
  lambda5 = 10,
  delta = 1,
  rho = 0.1,
  inclusion = 0.5,
  hom = c(3, 2)
)
```

## Arguments

- type:

  `"minnesota"`, `"ssvs"` or `"horseshoe"`.

- lambda1:

  overall tightness (also the SSVS slab scale).

- lambda2:

  cross-variable shrinkage factor (Minnesota).

- alpha:

  lag-decay exponent.

- lambda3:

  prior sd multiplier (times `s_i`) of deterministic coefficients.

- lambda5:

  prior sd multiplier of the contemporaneous coefficients \\a\_{ij}\\.

- delta:

  prior mean of the own first lag (1 random walk, 0 stationary).

- rho:

  SSVS spike-to-slab sd ratio.

- inclusion:

  SSVS prior inclusion probability.

- hom:

  prior `c(a0, b0)` for constant variances when `sv = FALSE`: IG(`a0`,
  `b0 s_i^2`).

## Value

a `prior_triangular` object.

## Details

The base prior variance of the coefficient on lag \\l\\ of variable
\\j\\ in equation \\i\\ is \\\frac{s_i^2}{s_j^2 l^{\alpha}}\\ (\\s\\ =
AR residual sd), multiplied by

- `minnesota`: \\\lambda_1^2\\ (own lags) or \\\lambda_1^2\lambda_2^2\\
  (other lags);

- `ssvs`: slab \\\lambda_1^2\\ or spike \\\rho^2\lambda_1^2\\ by a
  Bernoulli inclusion indicator;

- `horseshoe`: \\\tau_i^2\lambda\_{ij}^2\\ with half-Cauchy global/local
  scales (data-driven shrinkage).

## See also

Other Bayesian VAR:
[`autoplot.mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.md),
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md),
[`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md),
[`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.md),
[`tidy.bvarm()`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md),
[`volatility()`](https://kvasilopoulos.github.io/varm/reference/volatility.md)
