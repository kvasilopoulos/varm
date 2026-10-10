# Package index

## Reduced-form VAR

Estimation, lag selection and stability.

- [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md) :
  Estimate a reduced-form VAR by least squares
- [`var_select()`](https://kvasilopoulos.github.io/varm/reference/var_select.md)
  : Select the VAR lag order
- [`stability()`](https://kvasilopoulos.github.io/varm/reference/stability.md)
  : Eigenvalue moduli of the companion matrix
- [`lag_matrices()`](https://kvasilopoulos.github.io/varm/reference/lag_matrices.md)
  : Lag coefficient matrices

## Impulse responses, decompositions, forecasts

[`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md),
[`fevd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md),
[`hd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md) and
[`forecast()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
for VAR, BVAR and draws objects.

- [`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  [`fevd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  [`hd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  [`forecast()`](https://kvasilopoulos.github.io/varm/reference/dynamics.md)
  : Impulse responses, variance decompositions, historical
  decompositions, forecasts
- [`normalise_to()`](https://kvasilopoulos.github.io/varm/reference/normalise_to.md)
  : Normalise the impact of a shock

## Uncertainty

Bootstrap and Monte Carlo draws of the reduced form.

- [`draws()`](https://kvasilopoulos.github.io/varm/reference/draws.md) :
  Reduced-form parameter draws

## Structural identification

Recursive, long-run, sign, zero and sign, zero-pattern,
heteroskedasticity and user-supplied schemes.

- [`id_chol()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
  [`id_longrun()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
  [`id_given()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
  [`id_sign()`](https://kvasilopoulos.github.io/varm/reference/identification.md)
  : Structural identification schemes
- [`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.md)
  : Zero (short-run and long-run) restrictions in an A-B /
  Blanchard-Quah style
- [`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.md)
  : Sign restrictions combined with zero restrictions (Arias,
  Rubio-Ramirez & Waggoner 2018)
- [`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.md)
  : Identification through heteroskedasticity (Rigobon 2003), two
  regimes
- [`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.md)
  : Fry-Pagan (2011) median-target rotation

## Proxy-SVAR (external instruments)

Proxy identification, bootstrap, weak-instrument robust and Bayesian
inference.

- [`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.md)
  : External-instrument (proxy) SVAR identification

- [`proxy_fstat()`](https://kvasilopoulos.github.io/varm/reference/proxy_fstat.md)
  : First-stage strength of a proxy

- [`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.md)
  : Joint residual-proxy bootstrap for proxy SVARs

- [`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md)
  : Proxy-SVAR (external instrument) identification and bootstrap
  inference

- [`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.md)
  : Weak-instrument-robust confidence sets for proxy-SVAR impulse
  responses

- [`varm_proxy_bayes()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_bayes.md)
  :

  Bayesian proxy-SVAR on top of
  [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md)
  draws

- [`varm_proxy_arw()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_arw.md)
  [`id_proxy_bayes()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_arw.md)
  : Bayesian proxy-SVAR with exogeneity, relevance and sign restrictions
  (Arias-Rubio-Ramirez-Waggoner type)

- [`varm_proxy_gkr()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_gkr.md)
  : Robust-Bayes inference for set-identified proxy-SVARs (Giacomini,
  Kitagawa & Read)

## Bayesian VAR

Conjugate NIW (Minnesota, sum-of-coefficients, single-unit-root) and
triangular (shrinkage, stochastic volatility) samplers, with
diagnostics.

- [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.md) :
  Bayesian VAR
- [`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.md)
  : Minnesota / sum-of-coefficients / single-unit-root prior
- [`prior_triangular()`](https://kvasilopoulos.github.io/varm/reference/prior_triangular.md)
  : Triangular-sampler priors (large BVAR, shrinkage, stochastic
  volatility)
- [`volatility()`](https://kvasilopoulos.github.io/varm/reference/volatility.md)
  : Posterior volatility paths
- [`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md)
  [`tidy(`*`<mcmc_diagnostics>`*`)`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md)
  [`as_tibble(`*`<mcmc_diagnostics>`*`)`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md)
  [`glance(`*`<mcmc_diagnostics>`*`)`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.md)
  : MCMC convergence diagnostics

## Time-varying parameter VAR

TVP-VAR with stochastic volatility and its time-varying dynamics.

- [`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.md)
  : Time-varying parameter VAR with stochastic volatility
- [`tvp_prior()`](https://kvasilopoulos.github.io/varm/reference/tvp_prior.md)
  : Prior for the TVP-VAR (Primiceri training-sample calibration)
- [`irf(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md)
  [`fevd(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md)
  [`forecast(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.md)
  : Time-varying impulse responses, variance decompositions and
  forecasts of a TVP-VAR
- [`irf_path()`](https://kvasilopoulos.github.io/varm/reference/irf_path.md)
  : Response of one variable to one shock at a fixed horizon, across all
  dates
- [`tvp_volatility()`](https://kvasilopoulos.github.io/varm/reference/tvp_volatility.md)
  : Volatility paths of a TVP-VAR

## Local projections

Linear, IV, joint, state-dependent, smooth and panel local projections;
direct forecasts.

- [`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.md) : Local
  projections
- [`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.md)
  : Local projections with joint cross-horizon covariance
- [`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.md)
  : Sup-t simultaneous confidence bands for impulse responses
- [`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.md)
  : State-dependent local projections (Auerbach-Gorodnichenko)
- [`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.md)
  : Smooth local projections (Barnichon & Brownlees)
- [`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.md)
  : Panel local projections with fixed effects and clustered errors
- [`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.md)
  : Direct multi-step forecasts from local projections
- [`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.md)
  : Compare local-projection and VAR impulse responses

## Tidy output

[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`augment()`](https://generics.r-lib.org/reference/augment.html) methods
(tibbles with broom column names; `as_tibble()` is an alias of
[`tidy()`](https://generics.r-lib.org/reference/tidy.html)).

- [`tidy(`*`<varm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md)
  [`glance(`*`<varm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md)
  [`augment(`*`<varm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.md)
  : Tidy a VAR fit
- [`tidy(`*`<varm_irf>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_irf.md)
  [`tidy(`*`<varm_fevd>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_irf.md)
  [`as.data.frame(`*`<varm_dyn>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_irf.md)
  : Tidy impulse responses and variance decompositions
- [`tidy(`*`<varm_forecast>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_forecast.md)
  : Tidy forecasts
- [`tidy(`*`<varm_hd>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.varm_hd.md)
  : Tidy historical decomposition
- [`tidy(`*`<bvarm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md)
  [`as_tibble(`*`<bvarm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md)
  [`glance(`*`<bvarm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md)
  [`autoplot(`*`<bvarm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md)
  [`plot(`*`<bvarm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.md)
  : Tidy, glance at and plot a Bayesian VAR
- [`tidy(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
  [`as_tibble(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
  [`glance(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
  [`autoplot(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
  [`autoplot(`*`<tvpvar_path>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
  [`plot(`*`<tvpvar>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
  [`plot(`*`<tvpvar_path>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.md)
  : Tidy and plot a TVP-VAR
- [`tidy(`*`<lpm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`tidy(`*`<lp_state>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`tidy(`*`<lp_smooth>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`tidy(`*`<lp_forecast>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`tidy(`*`<lp_supt>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`as_tibble(`*`<lpm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`as_tibble(`*`<lp_state>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`as_tibble(`*`<lp_forecast>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`glance(`*`<lpm>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`glance(`*`<lp_panel>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`glance(`*`<lp_smooth>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  [`glance(`*`<lp_state>`*`)`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.md)
  : Tidy and glance at local-projection results

## Plots

[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods return ggplot2 objects;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
figure.
[`theme_varm()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)
and
[`varm_palette()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)
style them.

- [`autoplot(`*`<varm_irf>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.md)
  [`autoplot(`*`<varm_fevd>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.md)
  [`plot(`*`<varm_dyn>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.md)
  : Plot impulse responses and variance decompositions
- [`autoplot(`*`<varm_forecast>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_forecast.md)
  [`plot(`*`<varm_forecast>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_forecast.md)
  : Plot forecasts: history and fan chart
- [`autoplot(`*`<varm_hd>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_hd.md)
  [`plot(`*`<varm_hd>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_hd.md)
  : Plot historical decompositions
- [`autoplot(`*`<varm_draws>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_draws.md)
  [`plot(`*`<varm_draws>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_draws.md)
  : Plot parameter draws: trace and density
- [`autoplot(`*`<mcmc_diagnostics>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.md)
  [`plot(`*`<mcmc_diagnostics>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.md)
  : Plot MCMC diagnostics
- [`autoplot(`*`<lpm>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`autoplot(`*`<lp_state>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`autoplot(`*`<lp_smooth>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`autoplot(`*`<lp_supt>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`autoplot(`*`<lp_forecast>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`plot(`*`<lpm>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`plot(`*`<lp_state>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`plot(`*`<lp_supt>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  [`plot(`*`<lp_forecast>`*`)`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.md)
  : Plot local-projection results
- [`theme_varm()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)
  [`varm_palette()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.md)
  : varm plot theme and palette

## Data and re-exports

- [`yraw`](https://kvasilopoulos.github.io/varm/reference/yraw.md) :
  Quarterly US macroeconomic data
- [`reexports`](https://kvasilopoulos.github.io/varm/reference/reexports.md)
  [`tidy`](https://kvasilopoulos.github.io/varm/reference/reexports.md)
  [`glance`](https://kvasilopoulos.github.io/varm/reference/reexports.md)
  [`augment`](https://kvasilopoulos.github.io/varm/reference/reexports.md)
  [`autoplot`](https://kvasilopoulos.github.io/varm/reference/reexports.md)
  : Objects exported from other packages

## Package

- [`varm-package`](https://kvasilopoulos.github.io/varm/reference/varm-package.md)
  : varm: Fast Vector Autoregressions and Local Projections with an Rcpp
  Engine
