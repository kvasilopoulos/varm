# Structural identification schemes

Each constructor returns a small specification object consumed by
[dynamics](https://kvasilopoulos.github.io/varm/reference/dynamics.md),
and friends. The reduced-form innovations are \\u_t = P \varepsilon_t\\,
\\E\[\varepsilon_t\varepsilon_t'\]=I\\; an identification scheme picks
\\P\\ with \\PP' = \Sigma\\.

## Usage

``` r
id_chol(order = NULL)

id_longrun()

id_given(P)

id_sign(restrictions, horizon = 1L, nrot = 1L, maxtry = 1000L)
```

## Arguments

- order:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  variable ordering (first = most exogenous), e.g.
  `order = c(rate, infl, unemp)`; names, strings and positions all work.
  Defaults to the column order of the data.

- P:

  a `K x K` impact matrix.

- restrictions:

  a `K x K` matrix (rows: variables, columns: shocks) or a `K x K x H`
  array with entries `+1`, `-1` (sign) or `0` (unrestricted). Dimnames
  are used to label the shocks. Columns of all zeros are unidentified
  shocks.

- horizon:

  if `restrictions` is a matrix, the number of horizons (from impact)
  over which the restriction is imposed.

- nrot:

  number of accepted rotations to keep per reduced-form draw.

- maxtry:

  maximum rotations tried per reduced-form draw.

## Details

- `id_chol()` recursive (Cholesky) ordering, Sims (1980).

- `id_longrun()` Blanchard-Quah (1989): the cumulative long-run effect
  matrix is lower triangular.

- `id_sign()` sign restrictions on impulse responses (Uhlig 2005;
  Rubio-Ramirez, Waggoner & Zha 2010).

- `id_given()` user supplied impact matrix \\P\\.

## See also

Other identification:
[`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md),
[`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md),
[`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md),
[`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md),
[`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md),
[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)

## Examples

``` r
fit <- varm(yraw, p = 2)
id_chol(order = c(rate, infl, unemp))
#> <varm_id> recursive (Cholesky), order: c(rate, infl, unemp)
#> (first = most exogenous; the contemporaneous effect of later variables on earlier ones is zero)
irf(fit, id_chol(order = c(rate, infl, unemp)), horizon = 4, ci = "none")
#> <varm IRF> K = 3, horizon 0..4, id = chol

S <- matrix(c(1, 1, 1, 1, -1, 0, 1, 0, 0), 3, 3, byrow = TRUE)
id_sign(S)
#> <varm_id> sign restrictions (6 at impact; 1 horizon)
#>   + positive, - negative, · unrestricted; rows = variables, cols = shocks
#>   rotations: 1 per draw, max 1000 tries
#>    e1 e2 e3
#> v1  +  +  +
#> v2  +  -  ·
#> v3  +  ·  ·
```
