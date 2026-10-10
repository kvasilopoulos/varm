# Robust-Bayes inference for set-identified proxy-SVARs (Giacomini, Kitagawa & Read)

Prior-robust Bayesian inference on impulse responses in a proxy-SVAR
identified by exogeneity of `r` proxies plus (optional) sign
restrictions, following Giacomini, Kitagawa & Read (2022, *J.
Econometrics*) who build on Giacomini & Kitagawa (2021, *Econometrica*).
Write the SVAR in orthogonal-reduced form
\\y_t=Bx_t+\Sigma\_{tr}Q\varepsilon_t\\ and the first stage
\\m_t=Dy_t+Gx_t+\sum_l H_lm\_{t-l}+v_t\\. The reduced-form parameters
\\\phi=(B,\Sigma,D,G,H,\Upsilon)\\ are updated by the data; the rotation
\\Q\\ is not. The identifying restrictions are the exogeneity
restrictions \\D\Sigma\_{tr}q_j=0\\ for the first \\K-r\\ shocks (the
proxies are correlated with the **last `r` shocks**), the user's sign
restrictions \\S(\phi,Q)\ge0\\ and the sign normalisation
\\\mathrm{diag}(Q'\Sigma\_{tr}^{-1})\ge0\\. Instead of one prior for
\\Q\|\phi\\, GKR consider *all* priors supported on the set
\\\mathcal{Q}(\phi)\\ of admissible rotations; this gives, for each
\\\phi\\, bounds \\\[l(\phi),u(\phi)\]\\ of the response, and

- the **set of posterior means** \\\[E\_\phi l(\phi),\\E\_\phi
  u(\phi)\]\\, a consistent estimator of the identified set under strong
  proxies (GKR sec. 2.4), averaged over draws with a non-empty set;

- the **robust credible region** of level \\\alpha\\: the shortest
  interval that has posterior probability \\\ge\alpha\\ under every
  prior in the class, centred at \\\arg\min\_\eta \hat z\_\alpha(\eta)\\
  with radius \\\min\_\eta\hat z\_\alpha(\eta)\\, where \\\hat
  z\_\alpha(\eta)\\ is the \\\alpha\\-quantile over draws of
  \\\max\\\|\eta-l\|,\|\eta-u\|\\\\;

- the posterior probability that the identified set is non-empty
  (plausibility of the sign restrictions).

The bounds are approximated with GKR's **Algorithm 2** (Steps 2.1-2.3
and 3'): per reduced-form draw, up to `nq` rotations satisfying all
restrictions are drawn and the extremes taken. This underestimates the
identified set at each draw (converging as `nq` grows); GKR's
constrained-optimiser Algorithm 1 is **not** implemented. With `size`
set, Algorithm 3 is used (responses to a unit shock to `target`; bounds
may be unbounded if the impact on `target` can be near zero).
`lower`/`upper` of a point-identified shock (e.g. the proxy shock when
`r = 1`) coincide, and the set of posterior means then equals the
posterior mean of the point-identified response. The single-prior
posterior implied by the first accepted rotation of each draw (the
conventional Bayesian answer, GKR recommend reporting it too) is
returned as `draws`/`bands`.

## Usage

``` r
varm_proxy_gkr(
  fit,
  m,
  signs = NULL,
  shock = NULL,
  target = 1L,
  size = NULL,
  pm = 0L,
  horizon = 20L,
  level = c(0.68, 0.9),
  ndraw = NULL,
  nq = 100L,
  maxtry = 2000L,
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

- signs:

  sign restrictions on all shocks: `K x K` matrix or `K x K x H` array
  of `+1/-1/0` (rows variables, columns shocks; the last `r` columns are
  the proxy shocks), or `NULL`.

- shock:

  column (shock) whose responses are bounded; default `K` (the last
  proxy shock).

- target:

  variable(s) (name or position) whose impact response is normalised to
  `size`; length `r`.

- size:

  `NULL` for responses to a one-standard-deviation shock (GKR default),
  or the impact on `target` (single value) for responses to a unit shock
  (Algorithm 3).

- pm:

  number of lags of the proxy in the first-stage regression.

- horizon:

  maximum horizon (impact is 0).

- level:

  credibility levels of the robust credible regions.

- ndraw:

  number of posterior draws to use (default: all stored draws).

- nq:

  accepted rotations per reduced-form draw used for the bounds.

- maxtry:

  rotation draws attempted per reduced-form draw before declaring the
  set empty.

- sign_horizon:

  number of horizons (from impact) over which a *matrix* `signs` is
  imposed.

## Value

a `varm_proxy_gkr` object: `lower`, `upper` (set of posterior means,
`(H+1) x K`), `region` (`(H+1) x K x 2*length(level)`, robust credible
region limits ordered like `bands`), `prob_nonempty`, `lu` (list `l`,
`u` of `(H+1) x K x M` arrays over draws with non-empty set),
`draws`/`bands` (single-prior posterior, `shock` column only), and
metadata.

## Details

The reduced-form prior is the one in `fit` for \\(B,\Sigma)\\ (e.g.
Minnesota) and Jeffreys for the first stage (GKR sec. 4.1.1 use
independent Jeffreys priors; the NIW prior of
[`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
also fits their framework). Frequentist coverage is claimed by GKR only
for strong proxies (their sec. 3: with a weak proxy the Bernstein-von
Mises property fails and the set of posterior means need not be
consistent). Check
[`proxy_fstat()`](https://kvasilopoulos.github.io/varm/reference/proxy_fstat.md).
The conditioning number of \\D\Sigma\_{tr}\\ diagnostic of their sec.
4.1.6 is not computed.

## References

Giacomini, R., Kitagawa, T. and Read, M. (2022). Robust Bayesian
inference in proxy SVARs. *Journal of Econometrics* 228(1), 107-126.

Giacomini, R. and Kitagawa, T. (2021). Robust Bayesian inference for
set-identified models. *Econometrica* 89(4), 1519-1556.

## See also

[`varm_proxy_arw()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_arw.md),
[`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)
(frequentist weak-instrument-robust sets)

## Examples

``` r
# \donttest{
set.seed(2)
e <- matrix(rnorm(900), 300)
y <- matrix(0, 300, 3)
for (t in 2:300) y[t, ] <- 0.4 * y[t - 1, ] + e[t, ] %*% t(matrix(c(1, .3, .2, 0, 1, .4, .5, .2, 1), 3))
bf <- bvarm(y, p = 1, ndraw = 200, nburn = 100)
S <- matrix(0, 3, 3); S[1, 1] <- 1; S[2, 1] <- 1   # shock 1 raises variables 1 and 2 on impact
varm_proxy_gkr(bf, m = e[, 3] + rnorm(300, sd = 0.5), signs = S, shock = 1, horizon = 6, nq = 50)
#> <varm_proxy_gkr> robust-Bayes bounds, shock 1 (1 proxy; 1 s.d. shock), horizon 0..6
#> 200 reduced-form draws, nq = 50; P(identified set non-empty) = 1
#> set of posterior means [lower, upper] and 90% robust credible region:
#> y1:
#>     lower upper rcr_lo rcr_hi
#> h=0 0.285 1.027  0.185  1.097
#> h=4 0.014 0.043 -0.011  0.071
#> y2:
#>      lower upper rcr_lo rcr_hi
#> h=0  0.030 1.005  0.000  1.061
#> h=4 -0.004 0.061 -0.037  0.099
#> y3:
#>      lower upper rcr_lo rcr_hi
#> h=0  0.060 0.438  0.015  0.488
#> h=4 -0.007 0.011 -0.041  0.033
# }
```
