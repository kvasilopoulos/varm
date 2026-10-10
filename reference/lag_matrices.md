# Lag coefficient matrices

Lag coefficient matrices

## Usage

``` r
lag_matrices(object, ...)
```

## Arguments

- object:

  a `varm`/`bvarm` object or a coefficient matrix.

- ...:

  unused.

## Value

array `[K, K, p]` with `A[,,l]` the coefficient matrix on lag `l` (row =
equation, column = regressor variable).

## See also

Other estimation:
[`stability()`](https://kvasilopoulos.github.io/varm/reference/stability.md),
[`var_select()`](https://kvasilopoulos.github.io/varm/reference/var_select.md),
[`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)

## Examples

``` r
lag_matrices(varm(yraw, p = 2))[, , 1]
#>               infl      unemp        rate
#> infl   1.537775866 -0.2808590 -0.05713177
#> unemp -0.009517702  1.3296636 -0.02100610
#> rate   0.545670003 -0.6474938  0.77660355
```
