# Eigenvalue moduli of the companion matrix

A VAR is covariance stationary iff all moduli are \< 1.

## Usage

``` r
stability(object)
```

## Arguments

- object:

  a `varm` or `bvarm` object.

## Value

numeric vector of moduli (descending); for Bayesian objects, a matrix
`draws x Kp`.

## See also

Other estimation:
[`lag_matrices()`](https://kvasilopoulos.github.io/varm/reference/lag_matrices.md),
[`var_select()`](https://kvasilopoulos.github.io/varm/reference/var_select.md),
[`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)

## Examples

``` r
stability(varm(yraw, p = 2))
#>           [,1]
#> [1,] 0.9460450
#> [2,] 0.9460450
#> [3,] 0.7461457
#> [4,] 0.7461457
#> [5,] 0.4320615
#> [6,] 0.1394124
```
