
<!-- README.md is generated from README.Rmd. Please edit that file -->

# abvar

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/abvar)](https://cran.r-project.org/package=abvar)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/kvasilopoulos/abvar.svg?branch=master)](https://travis-ci.org/kvasilopoulos/abvar)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/kvasilopoulos/abvar?branch=master&svg=true)](https://ci.appveyor.com/project/kvasilopoulos/abvar)
[![Codecov test
coverage](https://codecov.io/gh/kvasilopoulos/abvar/branch/master/graph/badge.svg)](https://codecov.io/gh/kvasilopoulos/abvar?branch=master)
<!-- badges: end -->

## Introduction

The goal of {abvar} is the estimation and identification of Vector
Autoregressive (VAR) models. Still is a work in progess, any suggestion
for what might be of interest is much appreciated.

## Ideas on subroutines

Defers the computation until specific point

Naming: - `inst` - `gen` - `_`

## Status

  - [ ] VAR
      - [ ] VAR impulse response function (IRFs)
          - [ ] Identification
              - [ ] Reduce form
              - [ ] Cholesky
              - [ ] Long-run restrictions
              - [ ] Sign restrictions
                  - [ ] Uhlig
                  - [ ] Rubio Ramirez
                  - [ ] Fry Pagan
              - [ ] Proxy
                  - [ ] Mertens and Ravn; Stock and Watson
          - [ ] Confidence bands
              - [ ] Asymptotic
              - [ ] Monte Carlo
              - [ ] Bootstrap
              - [ ] Wild bootstrap
              - [ ] Block bootstrap
  - [ ] BVAR

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kvasilopoulos/abvar")
```

-----

Please note that the ‘abvar’ project is released with a [Contributor
Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.
