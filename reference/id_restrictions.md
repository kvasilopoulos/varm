# Zero (short-run and long-run) restrictions in an A-B / Blanchard-Quah style

General exactly- or under-identifying zero pattern on the impact matrix
\\P\\ (\\u_t=P\varepsilon_t\\) and on the long-run matrix
\\\Theta\_\infty = (I-\sum A_l)^{-1}P\\. Entries equal to `0` are
restricted to zero, `NA` (or anything else) are free. For a
just-identified pattern (after ordering shocks by number of zeros the
counts are \\K-1,K-2,\dots,0\\) the rotation is unique up to column
signs, found in closed form by sequential null-space projections (the
ARRW algorithm of
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md),
no numerical optimiser needed), so \\\Sigma=PP'\\ and every zero hold to
machine precision. Shock signs are normalised so that the impact
response of variable `j` to shock `j` (or the first free entry of column
`j`) is positive. Over-identified patterns are rejected (more zeros than
the counting allows); under-identified patterns give random draws from
the admissible set (see the missing-weights note in
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md)).

## Usage

``` r
id_restrictions(impact = NULL, longrun = NULL)
```

## Arguments

- impact, longrun:

  `K x K` matrices with `0` for zero restrictions and `NA` for free
  elements.

## Value

a `varm_id` object.

## Details

*Economic intuition.* The usual structural-VAR workhorse: Sims'
recursive scheme is `impact = lower triangular`, Blanchard-Quah is
`longrun = lower triangular`, and mixed schemes (e.g. Gali 1992,
Blanchard-Perotti with a long-run block) combine both. Pitfall: the
zeros are identifying assumptions that cannot be tested in the
just-identified case.

## See also

Other identification:
[`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md),
[`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md),
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md),
[`identification`](https://kvasilopoulos.github.io/varm/reference/identification.md),
[`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md),
[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)

## Examples

``` r
lower <- matrix(NA, 3, 3)
lower[upper.tri(lower)] <- 0
irf(varm(yraw, p = 2), id_restrictions(impact = lower), 3, ci = "none")
#> <varm IRF> K = 3, horizon 0..3, id = restrictions
#> sign-restriction acceptance rate: 1
```
