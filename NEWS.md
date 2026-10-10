# varm 0.1.0

First release. The package is experimental: the API may still change.

## 0.1.0 API notes

The API was refactored to tidy principles and ggplot2. Outputs now follow broom
column names; the object internals (`$estimate`, `$bands`, `$B`, `$draws`, ...)
are unchanged.

### Breaking or renamed outputs

* `as.data.frame()` on `irf()` / `fevd()` results (alias `tidy()`, `as_tibble()`):
  column `h` is now `horizon`; the wide `q0.xx` quantile columns are replaced by
  one row per cell and band `level` with `conf.low` and `conf.high`
  (`bands = FALSE` drops them). `response` and `shock` are factors in variable
  order.
* `lp_supt()` data frame / `tidy()`: columns are `response, horizon, estimate,
  std.error, band, level, conf.low, conf.high`, with `band` equal to `"sup-t"`
  or `"pointwise"`.
* `compare_lp_var()`: column `h` is now `horizon`.
* `mcmc_diagnostics()` tidy output: columns `param`, `mean`, `sd` are now `term`,
  `estimate`, `std.error` (plus `ess`, `rhat`, `geweke_z`).
* `irf_path()`: columns `median`, `lower`, `upper` are now `estimate`,
  `conf.low`, `conf.high` (with `date`).
* `plot()` methods return a `ggplot` (invisibly, after printing); base graphics
  are no longer used.

### New

* Generics re-exported from `generics` / `ggplot2`: `tidy()`, `glance()`,
  `augment()`, `autoplot()`. `tidy()` / `as_tibble()` methods for `varm`,
  IRF, FEVD, forecast, historical decomposition, `bvarm`, `tvpvar`, `lpm` and the
  `lp_*` results, `mcmc_diagnostics`; `glance()` for `varm`, `bvarm`, `tvpvar`,
  `lpm`, `lp_state`, `lp_smooth`, `lp_panel`, `mcmc_diagnostics`; `augment()` for
  `varm`.
* `autoplot()` methods for IRF, FEVD, forecasts, historical decomposition,
  draws, MCMC diagnostics, `bvarm`, `tvpvar` (coefficients, volatilities, IRF at
  dates) and every local-projection result.
* `theme_varm()`, `varm_palette()` (colour-blind safe) and `normalise_to()`
  (replaces the list form of `normalise`, which still works).
* tidyselect arguments for variable names (`lp(shock = , contemp = )`,
  `id_chol(order = )`, `autoplot(shocks = , responses = )`, ...), and a
  `vars =` argument of `varm()` for data frames.

### Fixes

* `bvarm()` works with `soc` / `sur` set to `NULL` in `prior_minnesota()`.
* `id_restrictions()`: the orientation of the solution is deterministic.
* `proxy_draws()`: revised block-length rule of the moving-block bootstrap.
* `draws(bias_correct = TRUE)` warns when combined with `method = "mc"` draws
  (it only applies to bootstrap draws).
* `irf.tvpvar()` rejects `id` and warns when a requested date is out of sample.
* `tvpvar()` warns when the prior degrees of freedom for `Q` do not exceed the number of coefficients + 1 (the prior mean of `Q` is then undefined, as with the default prior at `p = 4`).

## Reduced-form VAR

* `varm()`: OLS / Gaussian ML estimation of a VAR(p) with constant, linear
  trend and exogenous regressors; `var_select()` (AIC, HQ, SC, FPE);
  `stability()` (companion eigenvalue moduli); `lag_matrices()`; methods
  `print`, `summary`, `coef`, `resid`, `fitted`, `vcov`, `logLik`, `nobs`.
  Agrees with `vars` to about 1e-12 (coefficients, covariance, log-likelihood,
  information criteria, roots, Cholesky IRF, FEVD, point forecasts).
* `irf()`, `fevd()`, `hd()` and `forecast()` as S3 generics with methods for
  `varm`, `bvarm`, `varm_draws` (`irf`, `fevd`) and `tvpvar` (`irf`, `fevd`,
  `forecast`);
  unit-impact normalisation, cumulative responses, percentile bands at several
  levels, `autoplot()`, `plot()` and `tidy()` for IRF and FEVD
  objects.
* `draws()`: the common container (`varm_draws`) of reduced-form draws
  $(B, \Sigma)$ from the residual bootstrap, wild bootstrap (Rademacher or
  Mammen weights), Kilian bootstrap-after-bootstrap bias correction and the
  asymptotic normal / inverse-Wishart Monte Carlo.

## Structural identification

* `id_chol()`, `id_longrun()` (Blanchard-Quah), `id_sign()` (Haar rotations,
  Uhlig / Rubio-Ramirez-Waggoner-Zha), `id_given()`.
* `id_zerosign()`: sign plus zero restrictions (ARRW 2018 null-space sampler;
  the importance weights are not implemented, see the help page).
* `id_restrictions()`: short-run and long-run zero (A-B style) restrictions.
* `id_hetero()`: identification through heteroskedasticity, two regimes.
* `irf_fry_pagan()`: Fry-Pagan median-target rotation.
* Proxy-SVAR: `id_proxy()`, `proxy_draws()` (wild and moving-block bootstrap),
  `proxy_fstat()`, `varm_proxy()`, `varm_proxy_robust()` (weak-instrument robust
  confidence sets), `varm_proxy_bayes()` (limited-information Bayesian
  proxy-SVAR on `bvarm()` draws), `varm_proxy_arw()` (alias `id_proxy_bayes`,
  Arias-Rubio-Ramirez-Waggoner type posterior with exogeneity, relevance and
  sign restrictions) and `varm_proxy_gkr()` (Giacomini-Kitagawa-Read
  prior-robust bounds).

## Bayesian VAR

* `bvarm()` with `prior_minnesota()`: conjugate Normal-inverse-Wishart BVAR with
  Minnesota, sum-of-coefficients and single-unit-root priors; exact draws for
  fixed hyperparameters, optimiser plus adaptive random-walk Metropolis for the
  hyperparameters (Giannone-Lenza-Primiceri). The log posterior agrees with
  `BVAR::bvar` at matched prior scale to about 1e-13.
* `bvarm(prior = prior_triangular())`: triangular equation-by-equation Gibbs
  sampler with Minnesota, SSVS or horseshoe priors, optional stochastic
  volatility (Kim-Shephard-Chib mixture, forward-filtering backward-sampling);
  `volatility()`.
* `tvpvar()` with `tvp_prior()`, `irf_path()`, `tvp_volatility()` and
  `coef`, `summary`, `plot`, `irf`, `fevd`, `forecast` methods: TVP-VAR with
  stochastic volatility and time-varying contemporaneous matrix (Primiceri
  training-sample calibration).
* `mcmc_diagnostics()`: effective sample size, split R-hat, Geweke z-score for
  `bvarm` hyperparameters or user chains.

## Local projections

* `lp()`: local projections with Newey-West or HC1 standard errors, LP-IV
  (2SLS), lag augmentation, external shock series, cumulative responses.
  Matches `lm` + `sandwich` to about 1e-13.
* `lp_joint()` (joint cross-horizon covariance), `lp_supt()` (sup-t
  simultaneous bands), `lp_state()` (state-dependent LP; matches
  `lpirfs::lp_nl` to about 1e-12), `lp_smooth()` (smooth LP), `lp_panel()`
  (panel LP with fixed effects and clustered errors), `lp_forecast()` (direct
  forecasts) and `compare_lp_var()`.

## Data

* `yraw`: quarterly US inflation, unemployment and interest rate, 1953Q1-2006Q3
  (Koop and Korobilis).

## Engine

* All estimators, bootstraps, samplers and identification loops are compiled
  (Rcpp / RcppArmadillo, files in `src/`).

## Documentation and validation

* Vignettes: `getting-started`, `var-theory`, `identification`, `bayesian-var`,
  `tvp-var`, `local-projections`, `lp-advanced`, `choosing-a-method` and
  `validation` (reproducible validation report with live comparisons against
  `vars`, `BVAR`, `lpirfs`, `lm` + `sandwich` and quoted simulation studies).
* README rewritten with a feature matrix, benchmarks, validation table,
  limitations and a literature map; `_pkgdown.yml` lists every help topic once
  in grouped reference sections; `docs-src/methods-index.md` indexes every
  method with its definition and reference.
* Test suite of 65 tests / 728 expectations (`tests/testthat`) and heavier
  scripts in `validation/` (`bayes.R`, `ident-extra.R`, `lp-extra.R`,
  `proxy-coverage.R`, `tvp.R`; run the unit tests with
  `Rscript validation/run-all.R`).

## Known limitations

See the "Limitations" section of the README and the validation report:
under-coverage of draw-based IRF bands in short persistent or heteroskedastic
samples, the proxy-SVAR wild bootstrap at impact, sup-t bands, slightly low
stochastic-volatility calibration, slow mixing of the hierarchical-hyperparameter
sampler, unweighted zero+sign draws, and no missing-value handling.
