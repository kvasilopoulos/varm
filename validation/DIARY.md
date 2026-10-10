# varm validation diary — 2026-10-11

Scope: varm (kvasilopoulos/varm, HEAD d125cd9 + uncommitted changes) vs BEAR (8401f61), VAR-Toolbox (36e706a), BVAR_ (15bbc77).
Method caveat: no MATLAB/Octave available. "MATLAB side" = R ports of the .m source. BEAR ports were additionally checked against BEAR's shipped .mat results; VAR-Toolbox and BVAR_ ports only by internal consistency. Nothing was run in real MATLAB.
Details per toolbox: diary/baseline.md, diary/bear.md, diary/var-toolbox.md, diary/bvar.md. Scripts: validation/{bear,var-toolbox,bvar}/.

## Baseline
- Test suite: 652 PASS / 0 FAIL / 5 SKIP. `load_all` works on the installed build.
- Comparisons vs vars/BVAR/lpirfs/bsvars/svars are SKIPPED: those packages in rv/library don't load under Homebrew R (stale dylib paths). Unverified until reinstalled from source.
- Stale: tests/testthat/_problems/test-methods-26.R, testthat-problems.rds. `id_given` has no test.
- Working-tree src/RcppExports.cpp was regenerated (+232 lines) and `load_all` on the working tree did not compile at 00:24 (missing RcppArmadillo include) — check before committing.

## Agreement (varm = toolbox)
OLS B/Sigma, forecasts, Cholesky IRF, FEVD, HD, Blanchard-Quah, zero restrictions, sign-restriction distributions, NW marginal likelihood (no dummies), LP OLS/IV (matched settings), proxy-SVAR once conventions are matched.

## Suspected varm issues
| # | Issue | Source | Severity |
|---|-------|--------|----------|
| 1 | Marginal likelihood with soc/sur dummies wrong (-595.53 vs -630.38); feeds hierarchical sampler | BEAR | High |
| 2 | soc/sur dummies use first p regression rows, not presample (up to 0.14 in B). BVAR_ says this matches CRAN BVAR, so partly convention | BEAR, BVAR_ | Medium |
| 3 | `sur` dummy ignores exogenous regressors | BEAR | Low |
| 4 | Proxy zero-fill/uncentred moments: docs claim normalised IRFs unchanged; false unless instrument mean 0 (1.6-4.6% on GK). Demean on observed sample | VAR-Tbx, BVAR_ | Medium (doc/impl) |
| 5 | fevd() under id_proxy returns 100% for proxy shock (toolbox 17% at h=48) | VAR-Tbx | Misleading default |
| 6 | Bootstrap bands 6-7% wider (sqrt(Te/(Te-m)) rescale), undocumented; explosive draws not rejected | VAR-Tbx | Low (doc) |
| 7 | Normal-equation OLS loses ~6 digits on ill-conditioned level data (cpp_var_ols, cpp_lp) | VAR-Tbx, BVAR_ | Low |
| 8 | Default var=1e7 with soc=NULL shifts ML 0.3-8 log units vs flat constant | BVAR_ | Low |
| 9 | LP: default NW bandwidth h+1 (toolbox h); shock lags always added; normal not t quantiles | VAR-Tbx, BVAR_ | Low |

## Convention differences (not bugs)
AR variance divisor n-p vs n-2p-1; alpha = 2*lambda3; IW dof T-k vs T-k-K-1; forecast bands include parameter uncertainty; strict vs weak sign check; quantile type 7 vs 5; trend 1..T vs 1..Te; impact scale from full-sample vs instrument-sample Sigma; accept_rate counting.

## Bugs found in the toolboxes (for the record)
BVAR_: iresponse_longrun.m flips only Q(1,1); iresponse_heterosked.m uses singular values as covariance; global sign flip in sign check; zero+sign sampler non-uniform; "Jeffreys" label wrong; IC.m penalty too small; lambda=0 breaks ML.

## Not tested
No varm counterpart: narrative signs, ICA/higher-moment ID, FAVAR, mixed-frequency, panel, threshold, conditional forecasts, Pandemic/Long-run priors, BEAR Caldara-Herbst proxy. BEAR/BVAR_ SV/TVP models differ from varm's. IV bootstrap bands not compared numerically. Any real-MATLAB execution.

## Suggested next steps
1. Fix dummy-prior marginal likelihood and presample rows (#1-3); add tests against BEAR .mat values.
2. Fix proxy docs or demean instruments; revisit fevd under proxy.
3. Reinstall rv packages from source against Homebrew R to unskip comparison tests.
4. Decide on `brew install octave` to run toolbox code for real.

## Fix log — 2026-10-11 (dummy priors)
Edited src/bvar_niw.cpp, R/bvar.R (new args ybar/xbar), added src/varm_types.h (makes compileAttributes emit the RcppArmadillo include; it was being dropped, breaking compile), regenerated RcppExports.
- #2 soc/sur dummies now use the p presample observations (BEAR/Sims/GLP). BEAR T4a/b/c posterior tests now PASS (were FAIL, z up to 152). Old-convention T4a'/b'/c' checks removed from validation/bear/01_ols_nw.R.
- #3 `sur` dummy now puts presample mean/theta on all deterministic+exog columns (T4c PASS).
- #1 dummy-block ML now uses its posterior update. SumCoeff(0.1) logml: -631.3177 = BEAR formula (T11 PASS; was -595.5 with the stale rows). Residual: SumCoeff+InitialObs jointly differs by 0.045 (-612.869 vs -612.915), not yet explained — T11 FAIL; T11' "prior-mean" diagnostics are obsolete.
- Full suite (VARM_SRC=1): 846 PASS / 0 FAIL / 5 SKIP. Not committed. No new testthat tests added yet.
