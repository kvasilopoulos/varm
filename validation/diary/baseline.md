# Baseline audit - 2026-10-11

Environment: R 4.6.1 (Homebrew, arm64), rv library at `rv/library/4.6/arm64` (activated by `.Rprofile`). Nothing installed; no R sources edited; no commits.

## Package load
`pkgload::load_all(".")` - OK (Rcpp objects already compiled in `src/`). **PASS**

## Tests
`pkgload::load_all` + `testthat::test_dir("tests/testthat")` and `Rscript validation/run-all.R` (identical command): **[ FAIL 0 | WARN 0 | SKIP 5 | PASS 652 ]**, ~15 s. **PASS**

| file | expectations | fail | skip |
|---|---|---|---|
| test-bayes.R | 33 | 0 | 0 |
| test-bvar-conj.R | 30 | 0 | 1 (BVAR) |
| test-ident-extra.R | 339 | 0 | 0 |
| test-lp-extra.R | 28 | 0 | 1 (lpirfs) |
| test-lp.R | 38 | 0 | 0 |
| test-methods.R | 29 | 0 | 0 |
| test-psvar-bayes-full.R | 58 | 0 | 0 |
| test-psvar.R | 71 | 0 | 2 (`psvar` not installed) |
| test-tvp.R | 30 | 0 | 0 |
| test-var.R | (whole file skipped: vars) | 0 | 1 |

### Problems files
- `tests/testthat/_problems/test-methods-26.R` and `testthat-problems.rds` (both dated 2026-10-10 23:40): one recorded failure, "IRF object methods: print, as.data.frame, plot, autoplot" - `expect_output(print(id_chol()), "chol")` failed because the output was `<varm_id> recursive (Cho...`. **Stale**: test-methods.R passes now (29/29, 0 failures). The files are leftovers; safe to delete (not touched).

## Comparison against Suggests packages - **BLOCKED (not a varm failure)**
All of vars, BVAR, lpirfs, bsvars, svars, mvtnorm are present in `rv/library` but **cannot be loaded**; compiled deps link to CRAN-framework R paths that don't exist under Homebrew R:
- `lmtest.so` (needed by vars, lpirfs, svars): missing `/Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libgfortran.5.dylib`
- `mvtnorm.so` (BVAR), `stochvol.so` (bsvars): missing `.../libRblas.dylib`
Also `psvar` (not in Suggests) is absent.
Consequence: every cross-package comparison test is skipped (vars OLS/logLik/IC/stability in test-var.R; BVAR hierarchical NIW log-posterior; lpirfs::lp_nl; psvar). Those references are currently UNVERIFIED. Fix: `rv sync` / reinstall those packages from source against the Homebrew R (e.g. `rv` with source install, or symlink the missing dylibs) - not done here.
Standalone validation scripts `validation/{bayes,ident-extra,lp-extra,proxy-bayes-full,proxy-coverage,tvp}.R` were not run (run-all.R only wraps testthat).

## Exports vs tests
Method: word grep of the export name in `tests/testthat/*.R` (loose; may overcount generic names like `lp`, `draws`, `forecast`).
47 exports. Tested (by name): 46. **Not referenced in any test: `id_given`.**
Thinly covered (single file): compare_lp_var, id_hetero, id_proxy, id_proxy_bayes, id_restrictions, id_zerosign, irf_fry_pagan, irf_path, lag_matrices, lp_forecast, lp_joint, lp_panel, lp_smooth, lp_state, lp_supt, mcmc_diagnostics, proxy_draws, proxy_fstat, tvpvar family (tvp_prior, tvp_volatility; all in test-tvp.R), var_select, varm_proxy_arw/bayes/gkr/robust.
Multi-file: bvarm, draws, fevd, forecast, hd, id_chol, id_longrun, id_sign, irf, lp, prior_*, stability, varm, varm_proxy, volatility.
S3 methods (print/plot/summary/coef/etc.) are only exercised incidentally; test-methods.R is the sole dedicated file (4 tests).

## Verdict
- Load: PASS. Tests: PASS (652/0 fail). Cross-package validation: NOT VERIFIED (5 skips, env breakage). Coverage gap: `id_given`.
