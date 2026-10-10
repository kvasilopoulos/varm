# Bayesian proxy-SVAR on top of [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md) draws

A simple, limited-information Bayesian proxy-SVAR. For each posterior
draw \\(B_s, \Sigma_s)\\ of a
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
fit, with residuals \\u_t^{(s)}\\, the proxy is modelled by the
measurement regression \\m_t = \Lambda' u_t + e_t\\, \\e_t \sim N(0,
\Sigma_e)\\, with a flat (Jeffreys) prior on \\(\Lambda, \Sigma_e)\\. A
draw of \\(\Lambda_s, \Sigma\_{e,s})\\ from its (matrix-normal /
inverse-Wishart) conditional posterior gives \\\mathrm{Cov}(u, m) =
\Sigma_s\Lambda_s\\ and \\\mathrm{Var}(m) =
\Lambda_s'\Sigma_s\Lambda_s + \Sigma\_{e,s}\\, from which the impact
vector is obtained exactly as in
[`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md)
(also for `r > 1`) and the responses are \\\Theta_h(B_s)\\\gamma_s\\.

## Usage

``` r
varm_proxy_bayes(
  fit,
  m,
  target = 1L,
  size = 1,
  horizon = 20L,
  level = c(0.68, 0.9),
  ndraw = NULL
)
```

## Arguments

- fit:

  a [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
  object (class `bvarm`).

- m:

  proxy: numeric vector, or `T x r` matrix of `r` proxies, aligned with
  the data (the first `p` rows are dropped to match the residuals).

- target:

  variable(s) (name or position) whose impact response is normalised to
  `size`; length `r`.

- size:

  size of the shock (impact response of `target`); scalar or length `r`.

- horizon:

  maximum horizon (impact is 0).

- level:

  confidence levels; bands are equal-tailed percentiles.

- ndraw:

  number of posterior draws to use (default: all stored draws).

## Value

a `varm_proxy` object with `method = "bayes"` (same layout as
[`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md),
`draws` are posterior draws).

## Details

**Limitations (read before using).** This is *not* the joint-likelihood
Bayesian proxy-SVAR of Caldara & Herbst (2019) or Arias, Rubio-Ramirez &
Waggoner (2021): the proxy does not enter the likelihood of the VAR, so
it does not update \\(B,\Sigma)\\ and the reduced-form posterior is
whatever the
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
prior implies; there is no prior on the structural parameters or on the
instrument-shock relationship other than the flat prior on \\\Lambda\\;
stochastic-volatility fits (time-varying covariance) are rejected. With
a weak instrument the posterior of the normalised responses is correctly
wide (heavy-tailed) but, being Bayesian with a flat prior, its credible
sets do not have frequentist coverage guarantees (use
[`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)
for those). Credible bands are equal-tailed posterior percentiles;
`estimate` is the posterior median.

## References

Caldara, D. and Herbst, E. (2019). Monetary policy, real activity, and
credit spreads: evidence from Bayesian proxy SVARs. *AEJ:
Macroeconomics* 11(1), 157-192.

Arias, J. E., Rubio-Ramirez, J. F. and Waggoner, D. F. (2021). Inference
in Bayesian proxy-SVARs. *Journal of Econometrics* 225(1), 88-106.

Giacomini, R., Kitagawa, T. and Read, M. (2022). Robust Bayesian
inference in proxy SVARs. *Journal of Econometrics*
(doi:10.1016/j.jeconom.2021.02.003); not implemented.

## See also

[`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md),
[`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)

## Examples

``` r
# \donttest{
set.seed(1)
e <- matrix(rnorm(400), 200)
y <- matrix(0, 200, 2)
for (t in 2:200) y[t, ] <- 0.5 * y[t - 1, ] + c(e[t, 1] + 0.3 * e[t, 2], e[t, 2] + 0.5 * e[t, 1])
bf <- bvarm(y, p = 1)
varm_proxy_bayes(bf, m = e[, 1] + rnorm(200, sd = 0.5), horizon = 8, ndraw = 200)
#> <varm_proxy> shock to y1 (size 1), horizon 0..8, 200 posterior draws
# }
```
