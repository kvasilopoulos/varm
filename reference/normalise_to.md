# Normalise the impact of a shock

Rescales each structural shock so that the impact response of `variable`
to `shock` equals `size`. Pass the result to the `normalise` argument of
[`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md).

## Usage

``` r
normalise_to(shock, variable, size = 1)
```

## Arguments

- shock, variable:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  a single variable each: the shock being rescaled and the variable
  whose impact response is fixed.

- size:

  impact response of `variable` (default 1).

## Value

a `varm_normalise` specification.

## See also

Other dynamics:
[`irf.bvarm()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)

## Examples

``` r
normalise_to(rate, rate, size = 2)
#> $shock
#> <quosure>
#> expr: ^rate
#> env:  0x564a53777070
#> 
#> $variable
#> <quosure>
#> expr: ^rate
#> env:  0x564a53777070
#> 
#> $size
#> [1] 2
#> 
#> attr(,"class")
#> [1] "varm_normalise"
```
