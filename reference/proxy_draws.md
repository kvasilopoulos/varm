# Joint residual-proxy bootstrap for proxy SVARs

Draws `nboot` bootstrap samples of \\(B,\Sigma,\gamma)\\: `"wild"`
multiplies \\(\hat u_t, z_t)\\ by a common Rademacher weight
(Mertens-Ravn style; NOT recommended, since \\w_t^2=1\\ leaves
\\T^{-1}\sum u_t z_t\\ unchanged, so \\\gamma\\ only varies through
re-estimated VAR coefficients and bands are far too narrow, as shown by
Jentsch & Lunsford 2019 and in `validation/ident-extra.R`); `"mbb"`
(default) resamples moving blocks of the pairs \\(\hat u_t, z_t)\\
(block length `block`, default the Jentsch-Lunsford rule
\\\mathrm{round}(5.03\\T^{1/4})\\), preserving serial dependence in both
(Jentsch & Lunsford 2019). Each draw re-estimates the VAR on the
bootstrap data and recomputes \\\gamma^\* = T^{-1}\sum u^\*\_t
z^\*\_t\\. Percentile bands of the normalised IRF are only reliable with
a strong instrument.

## Usage

``` r
proxy_draws(
  fit,
  z,
  nboot = 1000L,
  type = c("mbb", "wild"),
  block = NULL,
  target = 1L,
  size = 1,
  shock = target
)
```

## Arguments

- fit:

  a [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)
  fit (used to compute the point `gamma`).

- z:

  proxy: numeric vector of length `nrow(fit$Y)` (aligned with the data;
  first `p` values dropped).

- nboot:

  number of draws.

- type:

  `"mbb"` (default, valid) or `"wild"` (invalid for gamma uncertainty,
  kept for comparison).

- block:

  block length for `"mbb"`.

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

a list with `draws` (a
[draws](https://kvasilopoulos.github.io/varm/reference/draws.md) object,
plus `gamma`), `id` (an
[`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md)
carrying one gamma per draw) and `fstat`. Use `irf(res$draws, res$id)`.

## See also

Other identification:
[`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md),
[`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md),
[`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md),
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md),
[`identification`](https://kvasilopoulos.github.io/varm/reference/identification.md),
[`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md)

## Examples

``` r
fit <- varm(yraw, p = 2)
z <- c(rep(0, 2), resid(fit)[, "rate"] + rnorm(nrow(resid(fit))))
pd <- proxy_draws(fit, z, nboot = 20, target = rate)
irf(pd$draws, pd$id, horizon = 4)
#> <varm IRF> K = 3, horizon 0..4, id = proxy, bands: 68%/90%
```
