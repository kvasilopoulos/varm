
<!-- README.md is generated from README.Rmd. Please edit that file -->

# varm

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/varm)](https://cran.r-project.org/package=varm)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/kvasilopoulos/varm/branch/master/graph/badge.svg)](https://codecov.io/gh/kvasilopoulos/varm?branch=master)
[![R-CMD-check](https://github.com/kvasilopoulos/varm/workflows/R-CMD-check/badge.svg)](https://github.com/kvasilopoulos/varm/actions)
<!-- badges: end -->

Fast frequentist and Bayesian vector autoregressions (VAR), structural
identification and local projections (LP) for R, with an
Rcpp/RcppArmadillo engine. Every sampler, bootstrap and identification
loop runs in compiled code.

The design idea is **one currency for uncertainty**: the bootstrap, the
asymptotic Monte Carlo and each Bayesian sampler return a `varm_draws`
object of reduced-form draws $(B, \Sigma)$. `irf()`, `fevd()`, `hd()`
and `forecast()`, and every identification scheme (`id_chol()`,
`id_longrun()`, `id_sign()`, `id_given()`), work on that object, so a
new estimator inherits all identification schemes and a new
identification scheme works with every estimator.

## Installation

``` r
# install.packages("pak")
pak::pak("kvasilopoulos/varm")
```

A C++ toolchain is required (the package compiles RcppArmadillo code).

## What is in the package

| Area | Methods | Entry points |
|----|----|----|
| Reduced-form VAR | OLS/ML, lag selection (AIC, HQ, SC, FPE), stability, summary | `varm()`, `var_select()`, `stability()` |
| Inference | residual / wild (Rademacher, Mammen) bootstrap, Kilian bias correction, asymptotic Monte Carlo | `draws()`, `irf(ci = )` |
| Dynamics | IRF, FEVD, historical decomposition, forecasts; unit-impact normalisation and cumulation | `irf()`, `fevd()`, `hd()`, `forecast()` |
| Identification | Cholesky, Blanchard-Quah long-run, sign restrictions (Haar rotations), user-supplied impact matrix | `id_chol()`, `id_longrun()`, `id_sign()`, `id_given()` |
| Bayesian VAR (conjugate) | Normal-inverse-Wishart with Minnesota, sum-of-coefficients and single-unit-root priors; closed-form marginal likelihood; hierarchical hyperparameters (GLP) | `bvarm()`, `prior_minnesota()` |
| Bayesian VAR (triangular) | Minnesota, SSVS and horseshoe priors; stochastic volatility (KSC mixture) | `bvarm(prior = prior_triangular())`, `volatility()` |
| Local projections | LP, LP-IV, lag augmentation, Newey-West / HC standard errors, cumulative responses | `lp()` |
| Proxy-SVAR | external-instrument identification with wild bootstrap | `varm_proxy()`, `proxy_fstat()` |

## Example

The `yraw` data set holds three quarterly US series (inflation,
unemployment, interest rate; Koop and Korobilis, 1953Q1-2006Q3).

``` r
library(varm)
y <- as.matrix(yraw); colnames(y) <- c("inflation", "unemployment", "rate")

fit <- varm(y, p = 2)
fit
#> <varm> VAR(2) with 3 variables, 213 effective observations
#> logLik = -338.416; max |eigenvalue| = 0.946
#> 
#>                 inflation unemployment    rate
#> inflation.l1       1.5378      -0.0095  0.5457
#> unemployment.l1   -0.2809       1.3297 -0.6475
#> rate.l1           -0.0571      -0.0210  0.7766
#> inflation.l2      -0.5418       0.0341 -0.4395
#> unemployment.l2    0.2336      -0.4290  0.6730
#> rate.l2            0.0580       0.0457  0.1232
#> const              0.2872       0.3654  0.0052
```

Responses to a policy-rate shock (recursive identification, rate ordered
last) with wild-bootstrap bands:

``` r
ir <- irf(fit, id = id_chol(), horizon = 16, ci = "wild", nboot = 300)
par(mfrow = c(1, 3), mar = c(3, 3, 2, 0.6), mgp = c(1.7, 0.5, 0))
for (v in colnames(y)) {
  b <- ir$bands[v, "rate", , ]; h <- 0:16
  plot(h, ir$estimate[v, "rate", ], type = "n", ylim = range(b), xlab = "horizon", ylab = "", main = paste("rate ->", v))
  polygon(c(h, rev(h)), c(b[, 1], rev(b[, 4])), border = NA, col = adjustcolor("steelblue", 0.2))
  polygon(c(h, rev(h)), c(b[, 2], rev(b[, 3])), border = NA, col = adjustcolor("steelblue", 0.3))
  abline(h = 0, col = "grey60"); lines(h, ir$estimate[v, "rate", ], lwd = 2, col = "steelblue4")
}
```

<img src="man/figures/README-irf-1.png" alt="" width="100%" />

The same call works for a Bayesian VAR (hierarchical Minnesota prior;
exact draws for the coefficients given the hyperparameters, Metropolis
for the hyperparameters):

``` r
bfit <- bvarm(y, p = 2, ndraw = 1000, nburn = 500)
bfit
#> <bvarm> niw-minnesota BVAR(2), K = 3, 1000 posterior draws
#> log marginal likelihood: -407.08849
#> hyperparameters (posterior mean [5%, 95%]); MH acceptance 0.512
#>          mean     5%    95%
#> lambda 0.6104 0.4069 0.8649
#> soc    0.4480 0.2036 1.2727
#> sur    0.6649 0.1953 1.5918
#> 
#> posterior mean of B:
#>                 inflation unemployment    rate
#> inflation.l1       1.5212       0.0224  0.5054
#> unemployment.l1   -0.2531       1.3314 -0.6012
#> rate.l1           -0.0451      -0.0200  0.7996
#> inflation.l2      -0.5267       0.0007 -0.4057
#> unemployment.l2    0.2159      -0.4082  0.6218
#> rate.l2            0.0458       0.0421  0.1055
#> const              0.2330       0.2514  0.0297
irb <- irf(bfit, id = id_chol(), horizon = 16)
round(irb$estimate["unemployment", "rate", c(1, 5, 9, 17)], 3)
#>     0     4     8    16 
#> 0.000 0.035 0.070 0.046
```

and a local projection of the same response:

``` r
lp(y, shock = "rate", horizon = 16, lags = 2, contemp = c("inflation", "unemployment"), lag_augment = TRUE)
#> <lpm> local projection of inflation, unemployment, rate on 'rate', horizons 0..16 (3 lags, HC)
#>    inflation unemployment    rate
#> 0     0.0000       0.0000  1.0000
#> 1    -0.0590      -0.0243  0.7763
#> 4    -0.0117       0.0456  0.6506
#> 8    -0.1030       0.2381  0.2098
#> 12   -0.5391       0.1760 -0.1398
#> 16   -0.6921       0.0634 -0.4354
```

## Learn more

- `vignette("var-theory")` - reduced-form VAR: definitions, IRF/FEVD/HD,
  bootstrap and Monte Carlo inference.
- `vignette("identification")` - Cholesky, long-run and sign
  restrictions; Fry-Pagan; normalisation.
- `vignette("bayesian-var")` - NIW posterior, Minnesota/unit-root
  priors, GLP hyperparameters, SSVS, horseshoe, stochastic volatility,
  diagnostics, comparison with `BVAR`.
- `vignette("local-projections")` - LP, LP-IV, lag augmentation, LP
  versus VAR bias-variance.
- `vignette("choosing-a-method")` - decision guide and a map of the
  VAR/LP literature and R packages (after
  [awesome-var](https://github.com/kvasilopoulos/awesome-var)).
- `docs-src/methods-index.md` - index of every method with its reference
  and where it is implemented.

## Related packages

`vars`, `svars`, `BVAR`, `bsvars`, `bsvarSIGNs`, `lpirfs`, `BigVAR`,
`bvartools`: see `vignette("choosing-a-method")` for how `varm` relates
to them. Results for the reduced-form VAR agree with `vars` to machine
precision, and the hierarchical BVAR agrees with `BVAR` (same prior
inputs) up to Monte Carlo error.

## Status

The package is experimental. See `NEWS.md` for the change log.
