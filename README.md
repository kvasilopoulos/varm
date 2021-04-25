
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

## Introduction

The goal of {varm} is the estimation and identification of Vector
Autoregressive (VAR) models. Still is a work in progess, any suggestion
for what might be of interest is much appreciated.

irf -&gt; array\[impulse, response, horizon, n\]

## Ideas on subroutines

Defers the computation until specific point

Naming: - `gen` - `_` - `impl`

## Status

-   [ ] VAR
    -   [x] VAR impulse response function (IRFs)
        -   [ ] Identification
            -   [ ] Reduce form
            -   [x] Cholesky
            -   [ ] Long-run restrictions
            -   [ ] Sign restrictions
                -   [ ] Uhlig
                -   [ ] Rubio Ramirez
                -   [ ] Fry Pagan
            -   [ ] Proxy
                -   [ ] Mertens and Ravn/Stock and Watson
        -   [ ] Confidence bands
            -   [ ] Asymptotic
            -   [ ] Monte Carlo
            -   [x] Bootstrap
            -   [ ] Wild bootstrap
            -   [ ] Block bootstrap
-   [ ] BVAR

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kvasilopoulos/varm")
```

------------------------------------------------------------------------

Please note that the ‘varm’ project is released with a [Contributor Code
of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.
