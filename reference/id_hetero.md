# Identification through heteroskedasticity (Rigobon 2003), two regimes

With two covariance regimes \\\Sigma_1=PP'\\ and \\\Sigma_2=P\Lambda
P'\\, \\\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_K)\\ with distinct
elements, \\P\\ is identified up to column order and sign: with
\\\Sigma_1=L_1L_1'\\, take the eigendecomposition
\\L_1^{-1}\Sigma_2L_1^{-T}=V\Lambda V'\\, \\P=L_1V\\. Shocks are ordered
by decreasing \\\lambda_j\\ (relative variance in regime 2), signs
chosen so \\\mathrm{diag}(P)\>0\\. The result is a *point* impact matrix
(\\\Sigma\\ from the full sample need not equal \\PP'\\; \\P\\ satisfies
it for the regime-1 covariance), returned as an
[`id_given()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
object with attribute `lambda`. It does not propagate sampling
uncertainty (all draws share one P).

## Usage

``` r
id_hetero(fit, regime)
```

## Arguments

- fit:

  a [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)
  fit.

- regime:

  logical or 0/1/2 vector of length `nrow(fit$Y)`, `TRUE`/`2` = regime
  2.

## Details

*Economic intuition.* Useful when volatility shifts (Great Moderation,
crisis episodes, announcement days) are plausibly shock-specific, with
no need for any economic restriction; the shocks are statistical and
must be labelled afterwards. Pitfalls: needs distinct variance ratios
(weak identification otherwise), regimes known/exogenous, and constant
\\P\\ across regimes.

## See also

Other identification:
[`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md),
[`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md),
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md),
[`identification`](https://kvasilopoulos.github.io/varm/reference/identification.md),
[`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md),
[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)

## Examples

``` r
fit <- varm(yraw, p = 1)
regime <- rep(0:1, length.out = nrow(yraw))
id_hetero(fit, regime)
#> <varm_id> user-supplied impact matrix P (response x shock)
#>        infl  unemp   rate
#> infl  0.192 -0.288 -0.146
#> unemp 0.140  0.273 -0.241
#> rate  0.476 -0.229  0.668
```
