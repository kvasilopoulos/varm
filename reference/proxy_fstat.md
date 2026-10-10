# First-stage strength of a proxy

Three F statistics for the relevance of the proxy for the `target`
shock: `F` from regressing the reduced-form residual of `target` on the
proxy, `F_robust` its heteroskedasticity-robust (HC1) version (needs
sandwich, else `NA`), and `F_hac`, the Montiel Olea-Stock-Watson (2021)
Wald statistic \\T\hat\Gamma_j^2/\hat W\_{\Gamma\_{jj}}\\ using the
delta-method HAC covariance that also accounts for the estimation of the
VAR coefficients. `F_hac` is the relevant statistic for
[`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md).

## Usage

``` r
proxy_fstat(fit, m, target = 1L, nwlags = NULL)
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

- nwlags:

  Newey-West truncation lag for `F_hac`; default
  `floor(4 * (T/100)^(2/9))`.

## Value

named numeric `c(F = , F_robust = , F_hac = )`. With several instruments
(`m` a matrix, `target` of the same length) a matrix with one row per
target; `F`/`F_robust` are then joint statistics over the instruments
and `F_hac` is `NA` (the MOSW sets are single-instrument).

## See also

[`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)
