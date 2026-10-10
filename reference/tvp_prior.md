# Prior for the TVP-VAR (Primiceri training-sample calibration)

Hyperparameters of the training-sample prior of Primiceri (2005). Let
\\\hat B,\hat\alpha,\hat\sigma\\ be OLS estimates on the first `tau`
observations (after lag loss), with covariance matrices \\V(\hat B),
V(\hat\alpha)\\. Then \$\$\beta_0\sim N(\hat B,\\ c\\V(\hat B)),\quad
\alpha_0\sim N(\hat\alpha,\\ c\\V(\hat\alpha)),\quad
\log\sigma\_{i0}\sim N(\log\hat\sigma_i, 10),\$\$ \$\$Q\sim
IW(k_Q^2\\\nu_Q\\V(\hat B),\nu_Q),\quad S_i\sim IW(k_S^2\\(1+\dim
S_i)\\V(\hat\alpha_i),1+\dim S_i),\quad W_i\sim IG(1, k_W^2).\$\$ The
`k` multipliers are the prior standard deviations of the *per-period*
drift relative to the OLS standard errors; they are the single most
influential choice in a TVP-VAR.

## Usage

``` r
tvp_prior(
  tau = 40L,
  kQ = 0.01,
  kS = 0.1,
  kW = 0.01,
  vscale = 4,
  df_Q = NULL,
  h_var = 10
)
```

## Arguments

- tau:

  training-sample length (effective observations, i.e. after the first
  `p`).

- kQ, kS, kW:

  prior scale multipliers for the innovations of the coefficients, of
  the contemporaneous (A) elements and of the log-volatilities.

- vscale:

  multiplier `c` on the OLS covariance in the prior for the initial
  states.

- df_Q:

  prior degrees of freedom of `Q` (default `tau`).

- h_var:

  prior variance of the initial log-variance (default 10).

## Value

a `tvp_prior` object for
[`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md).

## Details

**Economic intuition and pitfalls.** `kQ` encodes how fast the
transmission mechanism (lag coefficients) is believed to move. With
`kQ = 0.01` the prior says a coefficient drifts by about 1% of its OLS
standard error per quarter: over 170 quarters that is a standard
deviation of roughly 13% of the standard error times \\\sqrt{170}\\, so
substantial but gradual change is allowed. Raising `kQ` lets the data
attribute any misfit (outliers, omitted variables, mis-specified lags)
to coefficient drift, i.e. *spurious* time variation. The prior is not
neutral: with an inverse-Wishart prior whose degrees of freedom are tied
to a short training sample, the posterior of `Q` is dominated by the
prior in samples of typical macro length. Always report results for
several `kQ` values and compare with the constant-coefficient fit.

## See also

Other TVP-VAR:
[`irf_path()`](https://kvasilopoulos.github.io/varm/reference/irf_path.md),
[`tidy.tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md),
[`tvp_dynamics`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md),
[`tvp_volatility()`](https://kvasilopoulos.github.io/varm/reference/tvp_volatility.md),
[`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
