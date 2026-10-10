# Select the VAR lag order

Computes AIC, HQ (Hannan-Quinn), SC (Schwarz/BIC) and FPE for lags
`1..pmax` on a common estimation sample (the first `pmax` observations
are dropped for every `p`).

## Usage

``` r
var_select(Y, pmax = 8L, const = TRUE, trend = FALSE, exog = NULL, vars = NULL)
```

## Arguments

- Y:

  data frame (preferred), `ts` or numeric matrix (T x K) of endogenous
  variables. Non-numeric columns of a data frame are dropped; a `Date`
  column is kept as `$dates` and used by the plots.

- pmax:

  maximum lag.

- const:

  include a constant.

- trend:

  include a linear trend.

- exog:

  optional T x q matrix or data frame of exogenous regressors
  (contemporaneous, treated like deterministic terms; future values are
  then required for forecasting).

- vars:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  columns of `Y` that are the endogenous variables, e.g.
  `vars = c(infl, rate)`. Default: all numeric columns.

## Value

list with `selection` (the optimal lag per criterion) and the `criteria`
matrix.

## See also

Other estimation:
[`lag_matrices()`](https://kvasilopoulos.github.io/varm/reference/lag_matrices.md),
[`stability()`](https://kvasilopoulos.github.io/varm/reference/stability.md),
[`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)

## Examples

``` r
var_select(yraw, pmax = 4)$selection
#> AIC  HQ  SC FPE 
#>   3   2   2   3 
```
