# Changelog

## varm 0.1.0

First release. The package is experimental: the API may still change.

### Bayesian VAR fixes

- [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
  sum-of-coefficients and single-unit-root dummy priors now use the `p`
  presample observations (as in BEAR, Sims and GLP) instead of
  observations `p + 1` to `2p`.
- The single-unit-root dummy now covers exogenous and deterministic
  regressors, not only the constant.
- The marginal likelihood with dummy priors now updates the dummy
  block’s posterior; hierarchical `soc` / `sur` hyperparameter inference
  changes accordingly.

### 0.1.0 API notes

The API was refactored to tidy principles and ggplot2. Outputs now
follow broom column names; the object internals (`$estimate`, `$bands`,
`$B`, `$draws`, …) are unchanged.

#### Breaking or renamed outputs

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on
  [`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  /
  [`fevd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  results (alias
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  `as_tibble()`): column `h` is now `horizon`; the wide `q0.xx` quantile
  columns are replaced by one row per cell and band `level` with
  `conf.low` and `conf.high` (`bands = FALSE` drops them). `response`
  and `shock` are factors in variable order.
- [`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)
  data frame /
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html): columns
  are
  `response, horizon, estimate, std.error, band, level, conf.low, conf.high`,
  with `band` equal to `"sup-t"` or `"pointwise"`.
- [`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md):
  column `h` is now `horizon`.
- [`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md)
  tidy output: columns `param`, `mean`, `sd` are now `term`, `estimate`,
  `std.error` (plus `ess`, `rhat`, `geweke_z`).
- [`irf_path()`](https://kvasilopoulos.github.io/varm/reference/irf_path.md):
  columns `median`, `lower`, `upper` are now `estimate`, `conf.low`,
  `conf.high` (with `date`).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  return a `ggplot` (invisibly, after printing); base graphics are no
  longer used.

#### New

- Generics re-exported from `generics` / `ggplot2`:
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html),
  [`augment()`](https://generics.r-lib.org/reference/augment.html),
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html).
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  `as_tibble()` methods for `varm`, IRF, FEVD, forecast, historical
  decomposition, `bvarm`, `tvpvar`, `lpm` and the `lp_*` results,
  `mcmc_diagnostics`;
  [`glance()`](https://generics.r-lib.org/reference/glance.html) for
  `varm`, `bvarm`, `tvpvar`, `lpm`, `lp_state`, `lp_smooth`, `lp_panel`,
  `mcmc_diagnostics`;
  [`augment()`](https://generics.r-lib.org/reference/augment.html) for
  `varm`.
- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods for IRF, FEVD, forecasts, historical decomposition, draws,
  MCMC diagnostics, `bvarm`, `tvpvar` (coefficients, volatilities, IRF
  at dates) and every local-projection result.
- [`theme_varm()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md),
  [`varm_palette()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)
  (colour-blind safe) and
  [`normalise_to()`](https://kvasilopoulos.github.io/varm/reference/normalise_to.md)
  (replaces the list form of `normalise`, which still works).
- tidyselect arguments for variable names (`lp(shock = , contemp = )`,
  `id_chol(order = )`, `autoplot(shocks = , responses = )`, …), and a
  `vars =` argument of
  [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md) for
  data frames.

#### Fixes

- [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
  works with `soc` / `sur` set to `NULL` in
  [`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.md).
- [`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md):
  the orientation of the solution is deterministic.
- [`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md):
  revised block-length rule of the moving-block bootstrap.
- `draws(bias_correct = TRUE)` warns when combined with `method = "mc"`
  draws (it only applies to bootstrap draws).
- [`irf.tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md)
  rejects `id` and warns when a requested date is out of sample.
- [`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
  warns when the prior degrees of freedom for `Q` do not exceed the
  number of coefficients + 1 (the prior mean of `Q` is then undefined,
  as with the default prior at `p = 4`).

### Reduced-form VAR

- [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md):
  OLS / Gaussian ML estimation of a VAR(p) with constant, linear trend
  and exogenous regressors;
  [`var_select()`](https://kvasilopoulos.github.io/varm/reference/var_select.md)
  (AIC, HQ, SC, FPE);
  [`stability()`](https://kvasilopoulos.github.io/varm/reference/stability.md)
  (companion eigenvalue moduli);
  [`lag_matrices()`](https://kvasilopoulos.github.io/varm/reference/lag_matrices.md);
  methods `print`, `summary`, `coef`, `resid`, `fitted`, `vcov`,
  `logLik`, `nobs`. Agrees with `vars` to about 1e-12 (coefficients,
  covariance, log-likelihood, information criteria, roots, Cholesky IRF,
  FEVD, point forecasts).
- [`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md),
  [`fevd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md),
  [`hd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  and
  [`forecast()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  as S3 generics with methods for `varm`, `bvarm`, `varm_draws` (`irf`,
  `fevd`) and `tvpvar` (`irf`, `fevd`, `forecast`); unit-impact
  normalisation, cumulative responses, percentile bands at several
  levels,
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) for IRF and
  FEVD objects.
- [`draws()`](https://kvasilopoulos.github.io/varm/reference/draws.md):
  the common container (`varm_draws`) of reduced-form draws
  $`(B, \Sigma)`$ from the residual bootstrap, wild bootstrap
  (Rademacher or Mammen weights), Kilian bootstrap-after-bootstrap bias
  correction and the asymptotic normal / inverse-Wishart Monte Carlo.

### Structural identification

- [`id_chol()`](https://kvasilopoulos.github.io/varm/reference/identification.md),
  [`id_longrun()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
  (Blanchard-Quah),
  [`id_sign()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
  (Haar rotations, Uhlig / Rubio-Ramirez-Waggoner-Zha),
  [`id_given()`](https://kvasilopoulos.github.io/varm/reference/identification.md).
- [`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md):
  sign plus zero restrictions (ARRW 2018 null-space sampler; the
  importance weights are not implemented, see the help page).
- [`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md):
  short-run and long-run zero (A-B style) restrictions.
- [`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md):
  identification through heteroskedasticity, two regimes.
- [`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md):
  Fry-Pagan median-target rotation.
- Proxy-SVAR:
  [`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md),
  [`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)
  (wild and moving-block bootstrap),
  [`proxy_fstat()`](https://kvasilopoulos.github.io/varm/reference/proxy_fstat.md),
  [`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md),
  [`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)
  (weak-instrument robust confidence sets),
  [`varm_proxy_bayes()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_bayes.md)
  (limited-information Bayesian proxy-SVAR on
  [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
  draws),
  [`varm_proxy_arw()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_arw.md)
  (alias `id_proxy_bayes`, Arias-Rubio-Ramirez-Waggoner type posterior
  with exogeneity, relevance and sign restrictions) and
  [`varm_proxy_gkr()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_gkr.md)
  (Giacomini-Kitagawa-Read prior-robust bounds).

### Bayesian VAR

- [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
  with
  [`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.md):
  conjugate Normal-inverse-Wishart BVAR with Minnesota,
  sum-of-coefficients and single-unit-root priors; exact draws for fixed
  hyperparameters, optimiser plus adaptive random-walk Metropolis for
  the hyperparameters (Giannone-Lenza-Primiceri). The log posterior
  agrees with [`BVAR::bvar`](https://rdrr.io/pkg/BVAR/man/bvar.html) at
  matched prior scale to about 1e-13.
- `bvarm(prior = prior_triangular())`: triangular equation-by-equation
  Gibbs sampler with Minnesota, SSVS or horseshoe priors, optional
  stochastic volatility (Kim-Shephard-Chib mixture, forward-filtering
  backward-sampling);
  [`volatility()`](https://kvasilopoulos.github.io/varm/reference/volatility.md).
- [`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
  with
  [`tvp_prior()`](https://kvasilopoulos.github.io/varm/reference/tvp_prior.md),
  [`irf_path()`](https://kvasilopoulos.github.io/varm/reference/irf_path.md),
  [`tvp_volatility()`](https://kvasilopoulos.github.io/varm/reference/tvp_volatility.md)
  and `coef`, `summary`, `plot`, `irf`, `fevd`, `forecast` methods:
  TVP-VAR with stochastic volatility and time-varying contemporaneous
  matrix (Primiceri training-sample calibration).
- [`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md):
  effective sample size, split R-hat, Geweke z-score for `bvarm`
  hyperparameters or user chains.

### Local projections

- [`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md): local
  projections with Newey-West or HC1 standard errors, LP-IV (2SLS), lag
  augmentation, external shock series, cumulative responses. Matches
  `lm` + `sandwich` to about 1e-13.
- [`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md)
  (joint cross-horizon covariance),
  [`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)
  (sup-t simultaneous bands),
  [`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md)
  (state-dependent LP; matches
  [`lpirfs::lp_nl`](https://rdrr.io/pkg/lpirfs/man/lp_nl.html) to about
  1e-12),
  [`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md)
  (smooth LP),
  [`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md)
  (panel LP with fixed effects and clustered errors),
  [`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md)
  (direct forecasts) and
  [`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md).

### Data

- `yraw`: quarterly US inflation, unemployment and interest rate,
  1953Q1-2006Q3 (Koop and Korobilis).

### Engine

- All estimators, bootstraps, samplers and identification loops are
  compiled (Rcpp / RcppArmadillo, files in `src/`).

### Documentation and validation

- Vignettes: `getting-started`, `var-theory`, `identification`,
  `bayesian-var`, `tvp-var`, `local-projections`, `lp-advanced`,
  `choosing-a-method` and `validation` (reproducible validation report
  with live comparisons against `vars`, `BVAR`, `lpirfs`, `lm` +
  `sandwich` and quoted simulation studies).
- README rewritten with a feature matrix, benchmarks, validation table,
  limitations and a literature map; `_pkgdown.yml` lists every help
  topic once in grouped reference sections; `docs-src/methods-index.md`
  indexes every method with its definition and reference.
- Test suite of 65 tests / 728 expectations (`tests/testthat`) and
  heavier scripts in `validation/` (`bayes.R`, `ident-extra.R`,
  `lp-extra.R`, `proxy-coverage.R`, `tvp.R`; run the unit tests with
  `Rscript validation/run-all.R`).

### Known limitations

See the “Limitations” section of the README and the validation report:
under-coverage of draw-based IRF bands in short persistent or
heteroskedastic samples, the proxy-SVAR wild bootstrap at impact, sup-t
bands, slightly low stochastic-volatility calibration, slow mixing of
the hierarchical-hyperparameter sampler, unweighted zero+sign draws, and
no missing-value handling.
