# Sign restrictions combined with zero restrictions (Arias, Rubio-Ramirez & Waggoner 2018)

Let \\\Theta_h = \Phi_h L Q\\ be the structural IRF at horizon \\h\\,
with \\LL'=\Sigma\\ and \\Q\\ orthogonal, and \\\Theta\_\infty =
(I-\sum_l A_l)^{-1} L Q\\. A zero restriction on the response of
variable \\i\\ to shock \\j\\ at horizon \\h\\ reads \\e_i'\Phi_h L\\
q_j = 0\\, linear in the column \\q_j\\. Following ARRW Algorithm 2,
shocks are processed in decreasing order of their number of zero
restrictions \\z_j\\; column \\q_j\\ is drawn uniformly on the unit
sphere of the null space of \\\[Z_j L;\\ q_1';\dots;q\_{j-1}'\]\\, which
exists iff \\z\_{(j)} \le K-j\\. The draw is kept if all sign
restrictions (up to a per-shock sign flip) hold. Zeros hold exactly by
construction and \\PP'=\Sigma\\.

## Usage

``` r
id_zerosign(
  sign = NULL,
  zero = NULL,
  zero_lr = NULL,
  horizon = 1L,
  zero_horizon = 1L,
  nrot = 1L,
  maxtry = 1000L
)
```

## Arguments

- sign:

  `K x K` matrix or `K x K x H` array of `+1/-1/0` as in
  [`id_sign()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
  (`NULL` for none).

- zero:

  `K x K` matrix or `K x K x H` array with `1` where the response
  (variable x shock) is zero at that horizon (impact only when a matrix:
  `zero_horizon = 1`).

- zero_lr:

  `K x K` 0/1 matrix, `1` = long-run (cumulative) response is zero.

- horizon:

  horizon count over which a *matrix* `sign` is imposed.

- zero_horizon:

  horizon count for a *matrix* `zero` (default 1: impact only).

- nrot, maxtry:

  as in
  [`id_sign()`](https://kvasilopoulos.github.io/varm/reference/identification.md).

## Value

a `varm_id` object (usable in
[`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md) on
`varm_draws`/`varm`).

## Details

**Missing importance weights.** ARRW additionally reweight draws by a
Jacobian term (\\\|\det A_0\|^{-(2n+m+1)}\\ over \\\prod_i \|\det
V_i\|\\) so that they are uniform on the restricted set under the
conjugate prior. These weights are NOT implemented here: accepted draws
are *unweighted* samples from the null-space sampler, which is the exact
ARRW posterior only when no zero restrictions are present (then it
coincides with
[`id_sign()`](https://kvasilopoulos.github.io/varm/reference/identification.md))
or when the zero restrictions fully pin down the rotation up to sign
(then the weights are constant across rotations, though they still vary
over the reduced-form draws). Treat bands under partial zero-restriction
identification as approximate.

*Economic intuition.* Use it when theory gives a few exact zeros (e.g.
monetary policy has no contemporaneous effect on output and prices, or a
demand shock has no long-run effect on output) and signs for the rest.
Pitfalls: bands mix identification uncertainty and sampling uncertainty;
the implied prior over impulse responses is informative; zeros on many
shocks can make the problem exactly identified (sign information then
irrelevant).

## See also

[`id_sign()`](https://kvasilopoulos.github.io/varm/reference/identification.md),
[`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md),
[`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md)

Other identification:
[`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md),
[`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md),
[`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md),
[`identification`](https://kvasilopoulos.github.io/varm/reference/identification.md),
[`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md),
[`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)

## Examples

``` r
S <- matrix(0, 3, 3)
S[1, 1] <- 1
S[3, 3] <- 1
Z <- matrix(0, 3, 3)
Z[3, 2] <- 1
id <- id_zerosign(S, Z)
irf(draws(varm(yraw, p = 2), "mc", 20), id, horizon = 4)
#> <varm IRF> K = 3, horizon 0..4, id = zerosign, bands: 68%/90%
#> sign-restriction acceptance rate: 1
```
