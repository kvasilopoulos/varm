# Reduced-form parameter draws

Common currency for uncertainty: a `varm_draws` object holds `N` draws
of the coefficient matrix `B` (m x K x N) and innovation covariance
`Sigma` (K x K x N). Bootstrap, asymptotic Monte Carlo and every
Bayesian sampler in varm produce one, so identification, IRFs, FEVDs and
forecasts share one code path.

## Usage

``` r
draws(object, ...)

# S3 method for class 'varm'
draws(
  object,
  method = c("bootstrap", "wild", "wild_mammen", "mc"),
  n = 1000L,
  bias_correct = FALSE,
  ...
)
```

## Arguments

- object:

  a `varm` fit.

- ...:

  unused.

- method:

  `"bootstrap"` (recursive-design residual bootstrap), `"wild"`
  (Rademacher wild bootstrap, robust to conditional heteroskedasticity),
  `"wild_mammen"` or `"mc"` (asymptotic normal / inverse-Wishart Monte
  Carlo).

- n:

  number of draws.

- bias_correct:

  apply Kilian (1998) bootstrap-after-bootstrap bias correction.

## Examples

``` r
d <- draws(varm(yraw, p = 2), "mc", n = 50)
d
#> <varm_draws> 50 draws (mc); K = 3, p = 2
ggplot2::autoplot(d)
```
