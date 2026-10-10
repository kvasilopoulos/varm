# varm tidy design contract (v0.2)

Goal: tidy-tools style API, ggplot2 graphics, clean code. The numerical engine (src/*.cpp) and
the object *internals* (`$estimate`, `$bands`, `$B`, `$draws`, ...) stay as they are so every
existing validation test keeps passing; the tidy layer is added on top and the old
base-graphics code is removed.

## 1. Inputs
* Data: `data.frame`/tibble (preferred), `ts`, or matrix. Non-numeric columns (e.g. a `date`
  column) are dropped; a `Date`/`yearmon` column is kept as `$dates` for plots. Optional
  `vars = <tidyselect>` picks the endogenous columns (`varm(df, p = 2, vars = c(infl, rate))`).
* Every argument that names variables accepts **tidyselect**: bare names, strings, integers,
  `c()`, `all_of()`, `starts_with()`. One internal helper does it:
  `select_vars(quo, names, arg)` -> integer positions, with a clear `rlang::abort()` error.
  Applies to: `lp(shock=, contemp=, responses=, cumulative=)`, `id_chol(order=)`,
  `id_*` shock/variable arguments, `irf(cumulative=)`, `normalise_to(shock, variable, size)`
  (new, replaces the list argument; the old list form keeps working), `lp_state(state=)`, ...
  Capture with `rlang::enquo()` and resolve late (when the variable names are known).
* Errors/warnings/messages: `rlang::abort()/warn()/inform()` with a short message; no
  `stop(call. = FALSE)`.

## 2. Outputs: tidiers (broom conventions) -- all return tibbles
Generics re-exported from `generics`: `tidy()`, `glance()`, `augment()`; `tibble::as_tibble()`
methods are aliases of `tidy()` for the dynamic objects.
* `tidy.varm(x, conf.int = FALSE, conf.level = 0.95)`: `response, term, estimate, std.error,
  statistic, p.value` (+ `conf.low, conf.high`).
* `glance.varm`: one row: `K, p, nobs, logLik, AIC, BIC, HQ, max_root`.
* `augment.varm(x, data = NULL)`: the data plus `.fitted_<var>`, `.resid_<var>` (NA for the first p rows).
* `tidy.varm_irf` / `tidy.varm_fevd` (`as_tibble` alias): `response, shock, horizon, estimate`
  and, when bands exist, one row per (cell, level): `level, conf.low, conf.high`
  (`bands = FALSE` drops them). Factors: `response`/`shock` are factors in variable order.
* `tidy.varm_hd`: long: `time, response, shock, contribution` plus `baseline` rows
  (`shock == "baseline"`) and `actual` as its own column.
* `tidy.varm_forecast`: `horizon, variable, estimate (mean), level, conf.low, conf.high`.
* `tidy.bvarm`: posterior summaries of the coefficients: `response, term, estimate (mean),
  std.error (sd), conf.low, conf.high` (90% default, `conf.level =`); `glance.bvarm`:
  `K, p, nobs, ndraw, method, logml, accept` and hyper means as extra columns when present.
* `tidy.lpm` (and lp_state/lp_panel/lp_smooth/lp_joint results): `horizon, response, estimate,
  std.error, level, conf.low, conf.high` (+ `state` for lp_state).
* `tidy.tvpvar`: `date, variable/term, estimate, conf.low, conf.high`; volatilities likewise.
* `tidy.mcmc_diagnostics`, `glance` where meaningful.
Column names are lower-case snake/broom style; `horizon` (not `h`).

## 3. Graphics: ggplot2 only
* `autoplot()` (re-exported from ggplot2) methods for every result class return a `ggplot`:
  `varm_irf`, `varm_fevd` (facet_grid(response ~ shock), free y, nested ribbons by level),
  `varm_forecast` (history + fan chart), `varm_hd` (stacked columns/areas by shock + actual line),
  `lpm` and the lp_* results, `tvpvar` (coefficient/volatility paths; IRF at dates via
  `irf_path`), `bvarm` (`type = "coef"|"hyper"|"trace"`), `varm_draws`/`mcmc_diagnostics` (trace + density).
* `plot()` methods stay as thin wrappers that print the `autoplot()` and return it invisibly.
* One theme, `theme_varm(base_size = 11)` (minimal, light grid, bottom legend), one
  colour-blind-safe palette `varm_palette()`; no base graphics, remove `graphics`/`grDevices`
  from Imports. Axes labelled in words ("Horizon (periods)"), titles optional.

## 4. Code style (cleaner code)
* Run `air format .` (config in air.toml: 80 cols, 2-space) and keep `jarl` lint clean if available.
* Small single-purpose functions; no duplicated helpers; early returns; vectorised or compiled
  loops instead of nested R loops; no `;` chains; roxygen with `@family` tags and runnable `@examples`.
* Names: verbs for actions (`irf`, `fevd`, `forecast`, `draws`), `id_*` identification specs,
  `prior_*` priors, `lp_*` projections, `*_draws` posterior objects; nouns for results.
* Pipe-friendly: data/model is always the first argument; no hidden global state; deterministic
  given `set.seed`.
* Dependencies allowed: rlang, tibble, generics, ggplot2, tidyselect (+ Rcpp). No dplyr/tidyr in R/
  (keep it light); use base data manipulation internally and return tibbles.

## 5. Out of scope
src/*.cpp numerical code, the psvar*.R files (separately owned), vignettes (re-synced afterwards).
