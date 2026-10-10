# Fry-Pagan (2011) median-target rotation

Among the accepted draws of a set-identified
[`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
(sign / zero-sign), picks the single draw whose IRFs are closest, in
standardised squared distance over all responses, shocks and horizons,
to the pointwise median. It returns one internally consistent structural
model instead of the incoherent pointwise median.

## Usage

``` r
irf_fry_pagan(x)
```

## Arguments

- x:

  a `varm_irf` object with draws (from sign-type identification).

## Value

list with `index` (into `x$draws[,,,index]` and `x$draw`), `irf`
(`[K,K,H+1]`), `impact` (`K x K` impact matrix) and `varm_draw` (index
of the reduced-form draw it came from).

## Details

*Economic intuition.* Report this alongside the set bands when you need
one coherent shock series, historical decomposition or FEVD. Pitfall:
the median-target is itself a function of the prior over rotations and
is not "the" estimate.

## See also

Other identification:
[`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md),
[`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md),
[`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md),
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md),
[`identification`](https://kvasilopoulos.github.io/varm/reference/identification.md),
[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)

## Examples

``` r
S <- matrix(c(1, 1, 1, 1, -1, 0, 1, 0, 0), 3, 3, byrow = TRUE)
r <- irf(draws(varm(yraw, p = 2), "mc", 30), id_sign(S), horizon = 4)
irf_fry_pagan(r)$index
#> [1] 18
```
