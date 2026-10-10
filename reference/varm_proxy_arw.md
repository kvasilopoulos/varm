# Bayesian proxy-SVAR with exogeneity, relevance and sign restrictions (Arias-Rubio-Ramirez-Waggoner type)

Full-likelihood Bayesian proxy-SVAR in the spirit of Arias,
Rubio-Ramirez & Waggoner (2021). The `r` proxies obey the exogeneity
restrictions \\E\[m_t\varepsilon\_{jt}'\] = 0\\ for the first \\K-r\\
shocks and a relevance condition for the **last `r`** shocks. As in
their Proxy-SVAR (their eq. 3 with the block restrictions), the data
have the joint likelihood \\p(Y)\\p(M\mid Y)\\ with the proxy equation
\$\$m_t = D y_t + G x_t + \sum\_{l=1}^{p_m} H_l m\_{t-l} + v_t,\quad
v_t\sim N(0,\Upsilon),\$\$ which the lags of \\y_t\\ (and optionally of
\\m_t\\) are allowed to enter. The likelihood factorises, so a posterior
draw is: \\(B,\Sigma)\\ from the supplied
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md) fit
(any conjugate prior, e.g. Minnesota) and \\(D,G,H,\Upsilon)\\ from the
first-stage regression (Jeffreys prior), independently. The exogeneity
restrictions \\D\Sigma\_{tr}q_j=0,\\ j\le K-r\\, with
\\\Sigma=\Sigma\_{tr}\Sigma\_{tr}'\\ and \\Q\\ the rotation matrix, are
linear in the columns of \\Q\\; \\Q\\ is drawn column by column from the
unit sphere of the admissible subspace (ARW Algorithm 1 / Arias et al.
2018, using the orthogonal-reduced-form parameterisation), and sign
restrictions on the impulse responses of the proxy shocks are imposed by
accept/reject. With `r = 1` the proxy shock is **point-identified**
(\\q_K \propto \Sigma\_{tr}'D'\\) and signs only truncate the posterior;
with `r > 1` exogeneity alone leaves the proxy shocks identified only up
to an orthogonal \\r\times r\\ rotation, which the sign restrictions
(needed for identification) pin down. A reliability truncation
\\\lambda\_{\min}\\(\Upsilon+D\Sigma D')^{-1}D\Sigma D'\\\ge\gamma\\
(ARW 2021 section 2.3, after Caldara & Herbst 2019) is available through
`reliability`.

## Usage

``` r
varm_proxy_arw(
  fit,
  m,
  target = 1L,
  size = 1,
  signs = NULL,
  reliability = NULL,
  pm = 0L,
  horizon = 20L,
  level = c(0.68, 0.9),
  ndraw = NULL,
  nrot = 1L,
  maxtry = 1000L,
  sign_horizon = 1L
)

id_proxy_bayes(
  fit,
  m,
  target = 1L,
  size = 1,
  signs = NULL,
  reliability = NULL,
  pm = 0L,
  nrot = 1L,
  maxtry = 1000L,
  sign_horizon = 1L
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

  impact response of `target` to each proxy shock (recycled), or `NULL`
  for one-standard-deviation shocks (sign normalised so that \\A_0\\ has
  non-negative diagonal; state the sign through `signs` in that case).

- signs:

  sign restrictions on the responses of the proxy shocks: `K x r` matrix
  (or `K x r x H` array) of `+1`, `-1`, `0`; column `i` is the `i`-th
  proxy shock (shock `K - r + i`). `NULL` for none. The sign of a shock
  is otherwise normalised by `size` (impact on `target`).

- reliability:

  minimum reliability \\\gamma\in\[0,1)\\; draws below are discarded
  (truncated prior).

- pm:

  number of lags of the proxy in the first-stage regression.

- horizon:

  maximum horizon (impact is 0).

- level:

  confidence levels; bands are equal-tailed percentiles.

- ndraw:

  number of posterior draws to use (default: all stored draws).

- nrot:

  accepted rotations kept per reduced-form draw.

- maxtry:

  maximum rotation draws per reduced-form draw.

- sign_horizon:

  number of horizons (from impact) over which a *matrix* `signs` is
  imposed.

## Value

`varm_proxy_arw()`: a `varm_proxy` object (`method = "arw"`) with
`estimate` (posterior median), `draws`, `bands`, plus `accept_rate`,
`reliability` (posterior median of the reliability index, per proxy).
`id_proxy_bayes()`: a `varm_id` for
[`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md) on
a [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
fit (columns of non-proxy shocks are zero).

## Details

**What is and is not the ARW posterior.** The prior implied is the
conjugate prior on the reduced form times the (conditionally) uniform
prior on the admissible rotations, with importance weights set to one.
ARW's Algorithm 2 instead targets a normal-generalised-normal prior on
the *structural* parameters and re-weights draws by a volume element
(their Theorem 3, taken from Arias et al. 2018); that element is not
given in the paper and is **not implemented here**. The two posteriors
coincide when the proxy shocks are point-identified (`r = 1`, no sign
restrictions on non-proxy shocks) up to the reduced-form prior; with
`r > 1` they differ through the prior on the rotation. Non-proxy shocks
are not identified and their IRF columns are returned as zero. Compare
[`varm_proxy_bayes()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_bayes.md),
which does not allow lags of \\y_t\\ in the proxy equation and has no
sign or reliability restrictions, and
[`varm_proxy_gkr()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_gkr.md)
for prior-robust inference when `r > 1` or shocks other than the proxy
shocks are of interest.

Not implemented: Caldara & Herbst (2019)'s Metropolis-within-Gibbs
sampler with a structural measurement equation
\\m_t=\beta\varepsilon\_{1t}+\sigma\_\nu\nu_t\\ (where the proxy
likelihood depends on \\(B,\Sigma)\\ through \\u_t\\, so the proxy
updates the VAR coefficients) and priors on \\\beta,\sigma\_\nu\\.

## References

Arias, J. E., Rubio-Ramirez, J. F. and Waggoner, D. F. (2021). Inference
in Bayesian proxy-SVARs. *Journal of Econometrics* 225(1), 88-106.

Arias, J. E., Rubio-Ramirez, J. F. and Waggoner, D. F. (2018). Inference
based on SVARs identified with sign and zero restrictions.
*Econometrica* 86(2), 685-720.

Caldara, D. and Herbst, E. (2019). Monetary policy, real activity, and
credit spreads: evidence from Bayesian proxy SVARs. *AEJ:
Macroeconomics* 11(1), 157-192.

## See also

[`varm_proxy_gkr()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_gkr.md),
[`varm_proxy_bayes()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_bayes.md),
[`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md)

## Examples

``` r
# \donttest{
set.seed(1)
e <- matrix(rnorm(900), 300)
y <- matrix(0, 300, 3)
for (t in 2:300) y[t, ] <- 0.4 * y[t - 1, ] + e[t, ] %*% t(matrix(c(1, .3, .2, 0, 1, .4, .5, .2, 1), 3))
bf <- bvarm(y, p = 1, ndraw = 300, nburn = 100)
varm_proxy_arw(bf, m = e[, 3] + rnorm(300, sd = 0.5), target = 3, horizon = 8)
#> <varm_proxy> shock to y3 (size 1), horizon 0..8, 300 posterior draws
# }
```
