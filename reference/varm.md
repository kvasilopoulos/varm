# Estimate a reduced-form VAR by least squares

Equation-by-equation OLS (identical to Gaussian ML for the coefficients)
of \$\$y_t = c + \sum\_{l=1}^p A_l y\_{t-l} + \Gamma w_t + u_t,\quad u_t
\sim (0,\Sigma).\$\$ The estimation loop is compiled (RcppArmadillo).

## Usage

``` r
varm(Y, p = 1L, const = TRUE, trend = FALSE, exog = NULL, vars = NULL)
```

## Arguments

- Y:

  data frame (preferred), `ts` or numeric matrix (T x K) of endogenous
  variables. Non-numeric columns of a data frame are dropped; a `Date`
  column is kept as `$dates` and used by the plots.

- p:

  lag order.

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

an object of class `varm` with elements `B` (m x K coefficient matrix,
lags first then deterministics), `Sigma` (degrees-of-freedom corrected),
`U` residuals, `XtXinv`, `logLik` and metadata. See
[`tidy.varm()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md),
[`glance.varm()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md)
and
[`augment.varm()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md)
for tibble summaries.

## Details

See
[`vignette("var-theory")`](https://kvasilopoulos.github.io/varm/articles/var-theory.md)
for definitions and the economic reading of the estimates.

## See also

Other estimation:
[`lag_matrices()`](https://kvasilopoulos.github.io/varm/reference/lag_matrices.md),
[`stability()`](https://kvasilopoulos.github.io/varm/reference/stability.md),
[`var_select()`](https://kvasilopoulos.github.io/varm/reference/var_select.md)

## Examples

``` r
fit <- varm(yraw, p = 2)
fit
#> <varm> VAR(2) with 3 variables, 213 effective observations
#> logLik = -338.416; max |eigenvalue| = 0.946
#> 
#>             infl   unemp    rate
#> infl.l1   1.5378 -0.0095  0.5457
#> unemp.l1 -0.2809  1.3297 -0.6475
#> rate.l1  -0.0571 -0.0210  0.7766
#> infl.l2  -0.5418  0.0341 -0.4395
#> unemp.l2  0.2336 -0.4290  0.6730
#> rate.l2   0.0580  0.0457  0.1232
#> const     0.2872  0.3654  0.0052
summary(fit)
#> <varm> VAR(2), K = 3, T = 213
#> 
#> Equation: infl   R2 = 0.9833  adj.R2 = 0.9828
#>          Estimate Std. Error t value Pr(>|t|)    
#> infl.l1   1.53778    0.05965  25.780  < 2e-16 ***
#> unemp.l1 -0.28086    0.06218  -4.517 1.05e-05 ***
#> rate.l1  -0.05713    0.02955  -1.934 0.054527 .  
#> infl.l2  -0.54176    0.06131  -8.836 4.41e-16 ***
#> unemp.l2  0.23358    0.05920   3.946 0.000109 ***
#> rate.l2   0.05801    0.02906   1.996 0.047236 *  
#> const     0.28721    0.09550   3.008 0.002961 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Equation: unemp   R2 = 0.9487  adj.R2 = 0.9472
#>           Estimate Std. Error t value Pr(>|t|)    
#> infl.l1  -0.009518   0.064291  -0.148 0.882454    
#> unemp.l1  1.329664   0.067013  19.842  < 2e-16 ***
#> rate.l1  -0.021006   0.031845  -0.660 0.510224    
#> infl.l2   0.034078   0.066083   0.516 0.606629    
#> unemp.l2 -0.429014   0.063804  -6.724 1.71e-10 ***
#> rate.l2   0.045652   0.031324   1.457 0.146522    
#> const     0.365426   0.102927   3.550 0.000477 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Equation: rate   R2 = 0.9276  adj.R2 = 0.9255
#>           Estimate Std. Error t value Pr(>|t|)    
#> infl.l1   0.545670   0.148943   3.664 0.000316 ***
#> unemp.l1 -0.647494   0.155250  -4.171 4.47e-05 ***
#> rate.l1   0.776604   0.073776  10.526  < 2e-16 ***
#> infl.l2  -0.439517   0.153096  -2.871 0.004521 ** 
#> unemp.l2  0.673029   0.147816   4.553 9.03e-06 ***
#> rate.l2   0.123242   0.072568   1.698 0.090962 .  
#> const     0.005167   0.238454   0.022 0.982734    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual covariance:
#>            infl     unemp      rate
#> infl   0.094753 -0.003591  0.027003
#> unemp -0.003591  0.110071 -0.090048
#> rate   0.027003 -0.090048  0.590772
#> Residual correlation:
#>         infl  unemp   rate
#> infl   1.000 -0.035  0.114
#> unemp -0.035  1.000 -0.353
#> rate   0.114 -0.353  1.000
tidy(fit)
#> # A tibble: 21 × 6
#>    response term     estimate std.error statistic  p.value
#>    <chr>    <chr>       <dbl>     <dbl>     <dbl>    <dbl>
#>  1 infl     infl.l1   1.54       0.0596    25.8   2.13e-66
#>  2 infl     unemp.l1 -0.281      0.0622    -4.52  1.05e- 5
#>  3 infl     rate.l1  -0.0571     0.0295    -1.93  5.45e- 2
#>  4 infl     infl.l2  -0.542      0.0613    -8.84  4.41e-16
#>  5 infl     unemp.l2  0.234      0.0592     3.95  1.09e- 4
#>  6 infl     rate.l2   0.0580     0.0291     2.00  4.72e- 2
#>  7 infl     const     0.287      0.0955     3.01  2.96e- 3
#>  8 unemp    infl.l1  -0.00952    0.0643    -0.148 8.82e- 1
#>  9 unemp    unemp.l1  1.33       0.0670    19.8   1.09e-49
#> 10 unemp    rate.l1  -0.0210     0.0318    -0.660 5.10e- 1
#> # ℹ 11 more rows
glance(fit)
#> # A tibble: 1 × 8
#>       K     p  nobs logLik   AIC   BIC    HQ max_root
#>   <int> <int> <int>  <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1     3     2   213  -338.  731.  822.  768.    0.946

# choose the endogenous variables with tidyselect
varm(yraw, p = 1, vars = c(infl, rate))
#> <varm> VAR(1) with 2 variables, 214 effective observations
#> logLik = -356.646; max |eigenvalue| = 0.9798
#> 
#>            infl   rate
#> infl.l1  0.9971 0.0850
#> rate.l1 -0.0153 0.9044
#> const    0.0944 0.2006
```
