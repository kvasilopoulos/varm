# Proxy-SVAR (external instrument) identification and bootstrap inference

Mertens & Ravn (2013) / Stock & Watson (2012) proxy identification of
`r` structural shocks in a reduced-form
[`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md) fit,
using `r` external instruments. With one instrument the impact vector is
\\\gamma = T^{-1}\sum_t u_t m_t\\ (covariance of the reduced-form
residuals with the proxy), scaled so that the impact response of
`target` equals `size`. Responses are \\\Theta_h \gamma\\, reusing the
MA machinery of
[`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md).
With `r > 1` instruments the Mertens-Ravn (2013) `B11`/`B21` formulae
are used (see Details).

## Usage

``` r
varm_proxy(
  fit,
  m,
  target = 1L,
  size = 1,
  horizon = 20L,
  nboot = 500L,
  level = c(0.68, 0.9),
  method = c("mbb", "wild", "wild_normal"),
  block = NULL
)
```

## Arguments

- fit:

  a [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)
  object.

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

- nboot:

  number of bootstrap draws (`0` for point estimates only).

- level:

  confidence levels; bands are equal-tailed percentiles.

- method:

  `"mbb"` (default), `"wild"` or `"wild_normal"`; see Details.

- block:

  moving-block length for `method = "mbb"` (default Jentsch-Lunsford
  rule).

## Value

a `varm_proxy` object. For `r = 1`: `estimate` (`(H+1) x K`), `draws`
(`(H+1) x K x nboot`), `bands` (`(H+1) x K x 2*length(level)`). For
`r > 1` each gains a shock dimension after `K`. Also `gamma` (impact
covariance with the proxy), `fstat` (see
[`proxy_fstat()`](https://kvasilopoulos.github.io/varm/reference/proxy_fstat.md)),
`B1` (unit-variance impact matrix, `NULL` if not computable), `shocks`,
`reliability`, `method`, `block`.

## Details

**Inference (`method`).**

- `"wild"`:

  Mertens & Ravn (2013) / Gertler & Karadi (2015)-style recursive-design
  wild bootstrap: Rademacher weights multiply the residuals and the
  proxy jointly.

- `"wild_normal"`:

  same with standard-normal weights (the other variant studied by
  Jentsch & Lunsford).

- `"mbb"`:

  residual-based moving block bootstrap of Jentsch & Lunsford: blocks of
  length `block` are drawn jointly from the residuals and the proxy,
  each is recentred by the mean of the corresponding within-block
  position over all blocks, the VAR is re-simulated recursively and
  re-estimated. Block length default is the Jentsch-Lunsford rule \\\ell
  = 5.03\\T^{1/4}\\ (\\\ell = 20\\ at \\T = 250\\).

Jentsch & Lunsford (2019, 2022) show that the wild bootstrap is
asymptotically *invalid* for proxy-SVAR impulse responses, because
multiplying residuals and proxy by the same weight leaves the sample
covariance \\T^{-1}\sum m_t u_t'\\ (almost) fixed, so the sampling
variability of the impact vector is not reproduced, and that the moving
block bootstrap is valid. Mertens & Ravn (2019) dispute the practical
conclusions. The wild bootstrap remains available (`method = "wild"`) to
reproduce Mertens-Ravn / Gertler-Karadi bands. **The default changed
from `"wild"` to `"mbb"`** on the strength of that theory and of the
Monte Carlo in `validation/proxy-coverage.R` (wild-bootstrap 90%
coverage of the impact response of a non-target variable was about 14%
vs 87% for the block bootstrap). Neither bootstrap is robust to weak
instruments
([`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)
is). Intervals are equal-tailed percentiles, as in the papers.

**Several instruments.** `m` may be a `T x r` matrix and `target` a
vector of `r` variables: `r` instruments identify `r` shocks. As in the
Mertens-Ravn algorithm (and the SVAR-IV MATLAB code of
Miranda-Agrippino), the shocks are identified up to a Cholesky-type
restriction within the block of `target` variables, so the order of
`target` matters. The impact matrix is normalised so that column `i` has
impact `size[i]` on `target[i]`. Missing instrument values must be
zero-filled (this leaves the normalised responses unchanged); NA are
rejected.

The extra output `shocks` (realised unit-variance structural shocks,
\\\varepsilon_t = B_1'\Sigma_u^{-1}u_t\\) and `reliability` (squared
correlation of the instrument with the shock; squared canonical
correlations if `r > 1`) are sample objects computed from the point
estimates, valid when the instrument is relevant and exogenous; no
censoring adjustment is made.

## References

Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and
corporate income tax changes in the United States. *AER* 103(4),
1212-1247.

Jentsch, C. and Lunsford, K. G. (2019). The dynamic effects of personal
and corporate income tax changes in the United States: comment. *AER*
109(7), 2655-2678.

Jentsch, C. and Lunsford, K. G. (2022). Asymptotically valid bootstrap
inference for proxy SVARs. *Journal of Business & Economic Statistics*
(2022; working paper version: Cleveland Fed WP 16-19).

## See also

[`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md),
[`varm_proxy_bayes()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_bayes.md),
[`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md) (LP-IV)

## Examples

``` r
set.seed(1)
e <- matrix(rnorm(400), 200)
y <- matrix(0, 200, 2)
for (t in 2:200) y[t, ] <- 0.5 * y[t - 1, ] + c(e[t, 1] + 0.3 * e[t, 2], e[t, 2] + 0.5 * e[t, 1])
fit <- varm(y, p = 1)
varm_proxy(fit, m = e[, 1] + rnorm(200, sd = 0.5), horizon = 8, nboot = 100)
#> <varm_proxy> shock to y1 (size 1), horizon 0..8, 100 moving-block-bootstrap draws (block length 19)
#> first-stage F = 428.8; robust F = 441.7
#> instrument reliability: 0.766
varm_proxy(fit, m = e[, 1] + rnorm(200, sd = 0.5), horizon = 8, nboot = 100, method = "mbb")
#> <varm_proxy> shock to y1 (size 1), horizon 0..8, 100 moving-block-bootstrap draws (block length 19)
#> first-stage F = 520.8; robust F = 495.9
#> instrument reliability: 0.789
```
