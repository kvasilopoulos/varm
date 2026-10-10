# External-instrument (proxy) SVAR identification

`id_proxy()` identifies a single structural shock \\\varepsilon\_{1t}\\
from an external instrument \\z_t\\ (Mertens & Ravn 2013; Stock & Watson
2018; Gertler & Karadi 2015). With \\u_t = s\\\varepsilon\_{1t} +
\dots\\, the instrument satisfies relevance
\\E\[z_t\varepsilon\_{1t}\]=\phi\neq0\\ and exogeneity
\\E\[z_t\varepsilon\_{jt}\]=0,\\ j\>1\\, hence \$\$\gamma = E\[u_t z_t\]
= s\\\phi,\qquad s \propto \gamma .\$\$ The impact vector is identified
up to scale; it is normalised so that the impact response of variable
`target` equals `size`: \\s = \text{size}\cdot\gamma/\gamma\_{target}\\.
The response at horizon \\h\\ is \\\Phi_h s\\, \\\Phi_h\\ the
reduced-form MA coefficients. The estimate is the first stage of the
residual on the proxy: \\\hat\gamma = T^{-1}\sum_t \hat u_t z_t\\, and
relevance is judged by the first-stage F of the `target` residual on
\\z\\ (see
[`proxy_fstat()`](https://kvasilopoulos.github.io/varm/reference/proxy_fstat.md)).

## Usage

``` r
id_proxy(
  z = NULL,
  fit = NULL,
  gamma = NULL,
  target = 1L,
  size = 1,
  shock = target
)
```

## Arguments

- z:

  proxy: numeric vector of length `nrow(fit$Y)` (aligned with the data;
  first `p` values dropped).

- fit:

  a [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)
  fit (used to compute the point `gamma`).

- gamma:

  alternatively, a `K` vector or `K x N` matrix of reduced-form impact
  covariances (one column per reduced-form draw, as produced by
  [`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)).

- target:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variable whose impact response is normalised (name or position).

- size:

  impact response of `target` to the shock (e.g. 0.25 for a 25bp
  surprise).

- shock:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  column of the IRF array that holds the proxy shock (default `target`).

## Value

an object of class `varm_id`.

## Details

Only column `shock` of the structural impact matrix is identified: the
other columns of the returned IRF array are zero, so FEVD of the other
shocks is meaningless (the FEVD column of the proxy shock is correct
only if `size` is the shock standard deviation, which a proxy cannot
deliver).

Uncertainty needs the data, not only \\(B,\Sigma)\\: use
[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)
which bootstraps residuals and proxy jointly and carries one \\\gamma\\
per draw. Using `irf(fit, id_proxy(z, fit))` with `ci != "none"` holds
\\\gamma\\ fixed and therefore understates uncertainty.

*Economic intuition.* A macroeconomist uses a proxy when a credible
narrative/high-frequency measure of one shock is available
(monetary-policy surprises around FOMC announcements, tax narratives,
oil supply news) but a full recursive/sign scheme is not trusted.
Pitfalls: weak instruments (F below about 10 makes the normalised IRF
very wide and non-normal; use the Montiel Olea-Stock-Watson F_hac or
Anderson-Rubin bands, see
[`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)),
proxy contaminated by other shocks (exogeneity), and the sample used for
the proxy being shorter than the VAR (set missing proxy values to 0 only
if the instrument is truly 0 there).

## See also

[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md),
[`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md)

Other identification:
[`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md),
[`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md),
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md),
[`identification`](https://kvasilopoulos.github.io/varm/reference/identification.md),
[`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md),
[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)

## Examples

``` r
fit <- varm(yraw, p = 2)
z <- c(rep(0, 2), resid(fit)[, "rate"] + rnorm(nrow(resid(fit))))
id <- id_proxy(z, fit, target = rate)
irf(fit, id, horizon = 4, ci = "none")
#> <varm IRF> K = 3, horizon 0..4, id = proxy
```
