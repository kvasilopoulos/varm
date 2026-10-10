# varm vs ECB BEAR toolbox: validation diary (2026-10-11)

Scope: compare `varm` (installed build, `rv/library/4.6/arm64`, version 0.1.0, built 2026-10-11 00:08) with the ECB BEAR toolbox
(`/Users/VASILK04/Documents/research/external/BEAR-toolbox`, BEAR 5 core in `BEARX-Toolbox/bear/+bear/*.m`, BEAR 6 classes in
`BEARX-Toolbox/bear/bearing/+base/+estimator`). MATLAB/Octave are NOT available and were not installed. Method:

1. Read the MATLAB source of each overlapping method and replicate the deterministic parts in R (`validation/bear/common.R`).
2. Check the replicas against BEAR's own **shipped expected results**: `tests/bear5/results_ols.mat`, `results_test_data.mat`,
   `results_test_data_CH2019.mat` (MATLAB v5 files, parsed with the pure-Python reader `validation/bear/matread.py`; no numpy/scipy).
3. Compare varm to the replicas / stored BEAR numbers. Stochastic samplers are compared through Monte-Carlo z-scores or quantiles.

Data: BEAR tutorial `BEARX-tutorials/data/exampleData.csv` (== `data_.xlsx`): DOM_GDP, DOM_CPI, STN (+ Oil exogenous), 1974Q1-2014Q4,
p = 4, constant. BEAR convention: `startdate` includes the p presample rows, so the regression has T = 160 rows.

Scripts (run from the varm package root; each prints PASS/FAIL lines): `Rscript validation/bear/NN_*.R`
`01_ols_nw.R` (OLS, AR variances, Normal-Wishart posterior, dummies, log ML), `02_irf_fevd_hd.R`, `03_sign.R`, `04_stored_results.R`,
`05_zero.R`. `common.R` loads the **installed** varm (set `VARM_SRC=1` to use `pkgload::load_all` on the working tree).
Note: the working tree is shared with other sessions and currently does not recompile (see "Housekeeping" at the end).

Status legend: PASS = agreement within stated tolerance; FAIL = disagreement (classified below as bug or convention); NOT-TESTED = no counterpart / no reference.

## 1. Reduced-form OLS VAR

| id | what / how | numbers | status |
|---|---|---|---|
| T1 | `varm()` OLS B and Sigma vs BEAR `olsvar.m` replica, 3 endogenous + const + exogenous Oil, p=4 | max abs diff B 8e-14, Sigma 3e-16; Sigma uses T-k divisor in both | PASS |
| A1-A3 | varm vs **stored BEAR output** `results_ols.mat` (`Bhat`, `sigmahat`, `EPS`) | 5.5e-13, 4.4e-16, <1e-10 | PASS |
| C1-C2 | varm vs stored `results_test_data_CH2019.mat` (5 variables, p=12, 61 regressors, 1993m1-2007m6) | Bhat 2.9e-8 (ill-conditioned, 15-digit storage), Sigma ok | PASS |
| A4 | BEAR `olsforecast.m` 12-step point forecast (history ending 2013Q4, full-sample Bhat) vs `forecast.varm(ci="none")` on the same history | 1.7e-12 | PASS |
| A5 | BEAR forecast bands = point +- z * sqrt(cumulative sum Psi_j Sigma Psi_j') (parameters treated as known), replicated vs stored | 3e-13 | PASS (replica of BEAR) |
| A6 | varm `forecast(ci="mc")` bands vs BEAR band | varm 95% lower band h=1,6,12 (var 1): 0.605, -1.748, -2.138; BEAR 0.674, -1.429, -1.603 | convention difference: varm includes parameter uncertainty and (mc) draws, BEAR does not; not a bug |
| B8 | varm has no fixed-Sigma Minnesota (BEAR priors 11/12/13); see section 3 | - | NOT-TESTED |

## 2. Hyper-parameters / prior construction

BEAR NW prior (`nwprior.m`, prior 21): `B0` = ar on own first lag, `phi0` lag block var (1/arvar_i)(lambda1/j^lambda3)^2, constant/exogenous (lambda1*lambda4)^2,
`alpha0 = n+2`, `S0 = diag(arvar)`. varm `prior_minnesota` conjugate NIW: Omega_lag = lambda^2 / l^alpha / psi_i, IW(psi, K+2), prior mean `delta` on own lag, `var` for deterministic terms.
Exact mapping (used in all tests): `lambda = lambda1`, **`alpha = 2*lambda3`**, `psi = arvar`, `delta = ar`, `var = (lambda1*lambda4)^2`, `hyper = character(0)`.

| id | what / how | numbers | status |
|---|---|---|---|
| B1-B6 | my R replicas of `arloop.m`, `mprior.m`, `mpost.m` (prior 12) vs **stored** `arvar`, `beta0`, `omega0`, `sigma`, `betabar`, `omegabar` in `results_test_data.mat` | all < 1e-8 (typ. 1e-14) | PASS (validates my reading of BEAR; replicas then trusted) |
| T2 / B7 | AR(p)+const residual variance (the Minnesota/IW scale): varm `ar_resid_var` vs BEAR `arloop.m` (stored `arvar`) | ratio varm/BEAR = 0.96875 = (n-2p-1)/(n-p) exactly for all 3 series; BEAR divides by T-k (k=p+1), varm by n-p | FAIL = convention difference. varm's scale is ~3% smaller at p=4 and becomes larger for big p; not a bug but not BEAR-compatible. Worth documenting (or offering a df-corrected option) |
| conv-1 | lag-decay exponent: varm `alpha` is the exponent on the **variance** (lambda^2/l^alpha), BEAR `lambda3` is on the **sd** ((lambda1/j^lambda3)^2) | BEAR default lambda3=1 <=> varm alpha=2 | convention; the first run with alpha=1 failed by up to 339 sd, then PASS with alpha=2*lambda3 |
| conv-2 | defaults differ: varm lambda 0.2/alpha 2/delta 1/var 1e7/psi from n-p AR variances vs BEAR lambda1 0.1, lambda3 1, ar 0.8, lambda4 100 (var 100 on const) | - | convention |
| conv-3 | BEAR extras with no varm counterpart: lambda2 (cross-variable weight, Minnesota 11-13), lambda5 + block exogeneity, `Sigma="eye"` (S0 = I), exogenous prior means, per-variable ar vector via Excel | - | NOT-TESTED (feature gaps; GLP prior has no lambda2) |

## 3. Normal-Wishart posterior (fixed hyper-parameters)

Test: `bvarm(..., hyper = character(0), ndraw = 40000)` vs the analytic BEAR NW posterior (`nwpost.m` replica): E[B] = Bbar, E[Sigma] = Sbar/(alphabar-K-1), sd(B_ij) = sqrt(phibar_ii Sbar_jj/(alphabar-K-1)).

| id | case | result | status |
|---|---|---|---|
| T3 | NW, const + exogenous Oil, no dummies | E[B] max abs z 2.1 (diff 4e-4), E[Sigma] max z 1.2, sd ratio within 0.9% | PASS |
| T4a | + SumCoeff dummy (lambda6 = 0.1), BEAR dummies built from the p **presample** rows | E[B] z up to 136 (diff 0.14), sd ratio off 11% | FAIL (see bug 1) |
| T4a' | same but BEAR dummies built from the first p rows of the **regression** sample | z 2.5, sd within 1.3% | PASS (this is what varm does) |
| T4b / T4b' | + InitialObs (lambda7 = 1e-3), no exog; presample rows vs regression rows | FAIL (z 142, diff 0.05) / PASS (z 1.5) | as above |
| T4c / T4c' | + InitialObs with exogenous Oil | FAIL even with regression rows: z 93, diff 0.07 | FAIL (see bug 2) |
| T11 | log marginal likelihood: BEAR `nwmlik.m`/`nwmlikdata.m` formula vs `bvarm(...)$logml`, no dummies | -630.4095731 both | PASS |
| T11 | same with SumCoeff dummy | varm -595.5317, BEAR -630.3788 | FAIL (see bug 3) |
| T11'' | independent general NW formula for p(Y \| dummy-updated prior) | -630.3787685 = BEAR to 1e-6 | confirms BEAR value is the right one |
| T11' | BEAR formula but with the dummy block evaluated at the prior mean (no posterior update) | matches varm to 9e-11 (SumCoeff only) | explains bug 3 |

Not tested: BEAR Minnesota priors 11/12/13 (Sigma fixed), Independent-NW (prior 31/32 Gibbs), Normal-diffuse (41), flat (51), BEAR "hogs" optimiser and Minnesota-dummy / long-run dummies (no varm counterpart; varm hyper-optimisation is GLP-style and uses the same ML formula tested in T11, so it agrees with BEAR for fixed hyper and no dummies). NOT-TESTED.

## 4. IRF / FEVD / HD (Cholesky)

BEAR `irfsim.m` recursion (shock vector at lag p, `X = [lags, 0]` so exogenous terms are zero), `D = chol(Sigma,'lower')`, `fevd.m` (cumulated squared orthogonalised IRFs / total), `hdecomp.m` (sum_l ortirf(l) * shock(t-l+1), shocks = D^-1 u, baseline = Y - sum).

| id | what | numbers | status |
|---|---|---|---|
| T5 | orthogonalised IRF, 13 horizons, 3x3, point estimate | 2.6e-13 | PASS |
| T6 | Cholesky with reordering (STN first) | 3.0e-13 (replica of ordering semantics is my construction, BEAR `Cholesky.m` reorders then back-orders columns) | PASS |
| T7 | FEVD 13 horizons, rows sum to 1 | 1.1e-13 | PASS |
| T8 | HD contributions (160 x 3 x 3) and baseline | 3.3e-12 / 3.9e-12 | PASS |
| - | BEAR `hdecomp_inc_exo` splits baseline into deterministic and exogenous parts; varm returns one baseline | - | feature difference |
| - | BEAR IRFt=3 (triangular factorisation F Lambda F'), long-run IRF from `irfsim` with horizon 1000 | varm `id_longrun`/`id_chol` exist but BEAR long-run path not replicated | NOT-TESTED |

## 5. Sign / zero / relative-magnitude restrictions

| id | what / how | numbers | status |
|---|---|---|---|
| T9 | BEAR `irfres.m`+`checksignres.m` sequential column sampler (replica) vs `id_sign(nrot=3000)` at fixed OLS (B, Sigma); 3 shocks, 4 impact-horizons | acceptance rate 0.0119 (varm) vs 0.0120 (replica); max z of mean IRF difference 2.3; max quantile diff 0.12 sd | PASS |
| T9b | every accepted varm rotation satisfies all signs | 0 violations / 3000 | PASS |
| T10 | **stored BEAR draws** `D_record` (1001 sign-restricted impact matrices, `results_ols.mat`; table Demand(+,+,+) Supply(GDP+,CPI-) Money(GDP+,CPI+,STN-) + relative-magnitude GDP_dem >= GDP_sup, parameters from `olsirft4.m`: Sigma^-1 ~ W(Sigma^-1/T, T), beta ~ N) vs varm `ci="mc"`, 60000 draws, `id_sign` + post-hoc relmagn filter (31438 kept) | medians match within 0.10 sd (e.g. GDP<-Dem 0.638 vs 0.640, STN<-Mon -0.357 vs -0.362), 5/95% quantiles within 0.017 | PASS |
| T13 | BEAR `find_medianmodel.m` (standardised) vs varm Fry-Pagan on 3000 draws | same index (612) | PASS |
| T12 | recursive impact zeros + positive diagonal reproduce Cholesky exactly (`id_zerosign`) | maxdiff 0 | PASS |
| T12b | varm zero restrictions hold exactly for permuted shock order | PASS | PASS |
| - | BEAR `qzerores.m` draws columns in user order; with the unrestricted shock first the later columns have no admissible null space (replicated); varm/ARRW sorts shocks by number of zeros | - | note (BEAR limitation, not varm) |
| - | BEAR magnitude (`Mcell`) and relative-magnitude restrictions, FEVD restrictions, correlation-with-instrument restrictions | no varm counterpart (varm only filtered post hoc in T10) | NOT-TESTED |
| - | strict (>) vs weak (>=) sign check: varm `check_col` rejects zeros, BEAR `checksignres.m` accepts `>=0` | measure-zero | convention |
| - | importance weights for zero+sign restrictions (ARRW Jacobian) absent in varm; BEAR does not weight either | - | same as BEAR |

## 6. Proxy / IV

BEAR IRFt=6 (frequentist Gertler-Karadi, `irfiv_ols_for_bootstrap_GK.m`): 2SLS u_h on fitted u_1|[1,z], `b_1=1`, `D1 = b / ||chol(Sigma_IV)^-1 b||`, `Sigma_IV = EPS'EPS/(T-k)`. BEAR IRFt=5 (`irfIV_MH.m`, Caldara-Herbst Bayesian proxy SVAR).

| id | what / how | numbers | status |
|---|---|---|---|
| C3 | CH2019 data (`data_CH2019.xlsx`, MHF instrument, p=12, 5 variables, 162 obs): BEAR GK 2SLS impact (replica) x Psi_h vs `varm_proxy(size=1)` 49 horizons | 6e-6 | PASS |
| C4 | unit-variance impact vector: varm `B1` / BEAR `D1` | ratio 0.7896 for all 5 elements = sqrt((T-k)/T) exactly (BEAR Sigma/(T-k), varm Sigma/T) | PASS as convention |
| C5 | BEAR IRFt=5 stored posterior (`D_record`, 1800 draws, `results_test_data_CH2019.mat`) vs varm | varm documents that the Caldara-Herbst Metropolis-within-Gibbs sampler is NOT implemented; `varm_proxy_bayes/arw/gkr` are different posteriors. Info only: stored posterior-median D[,1] = (0.136, -0.015, -0.003, 0.011, 0.083), frequentist varm B1 = (0.116, -0.059, -0.004, -0.140, 0.054) (same sign pattern for 3 of 5; the Bayesian one has the prior and the proxy updating the VAR) | NOT-TESTED |
| - | BEAR F-test of instrument relevance (`irfiv_ols.m`) vs `proxy_fstat` | formula identical (regress residual on [1, z], classical F); not run numerically | NOT-TESTED |

## 7. Stochastic volatility, time-varying parameters, other BEAR models

NOT-TESTED, and not comparable by design:

* BEAR SV models: `CogleySargentSV` (stvol1, per-variable AR(1) log-vol, MH), `RandomInertiaSV` (stvol2), `CarrieroSV` (stvol3: **one common** scalar log volatility, MH). varm SV = triangular equation-by-equation sampler with K independent random-walk log-vols (Kim-Shephard-Chib mixture). Different models, no stored BEAR output, stochastic posteriors.
* BEAR `BetaTV`/`GeneralTV` (stvol/tvbvar 1-2: random-walk betas with Q ~ IG(chi=200, psi=0.01), tau=10000) vs varm `tvpvar` (Primiceri 2005 training-sample prior, kQ/kS/kW). Different priors.
* BEAR panel, FAVAR, mixed-frequency, threshold, mean-adjusted, conditional forecasts, dummy-observation Minnesota/long-run priors: no varm counterpart.
* varm-only features (local projections, LP-IV, bootstrap/wild/mbb bands, robust proxy inference) have no BEAR counterpart.

## Discrepancies and suspected varm bugs

1. **Dummy-observation prior built from the wrong rows (suspected bug)** `src/bvar_niw.cpp`: `make_spec` sets `ybar = mean(Y.rows(0, p-1))` where `Y` is the regression sample (`ols$Y`, i.e. observations p+1..2p) rather than the p initial/presample observations (BEAR `dummies.SumCoeff/InitialObs`: `longY(1:order,:)`; Sims / GLP use the first p observations y_1..y_p). Evidence: T4a/T4b FAIL against BEAR; replicating BEAR but taking rows p+1..2p (T4a'/T4b') PASSES. Affects `soc` and `sur` posteriors (up to 0.14 in B).
2. **`sur` dummy ignores exogenous regressors (suspected bug / limitation)**: the initial-observation dummy puts 0 on exogenous columns (only the constant gets 1/theta); BEAR uses mean of the initial exogenous values / lambda7 (T4c, T4c' FAIL, z 93). Only matters when `exog` is used with `sur`.
3. **Marginal likelihood with dummy priors is not p(Y | dummies) (suspected bug)**: the correction term for the dummy block is evaluated at the prior mean without the posterior update (`ev(Yd, Xd, est = false)` gives `bh = b`, `pen = 0`). The NW marginal likelihood of the dummy sample needs the updated posterior (BEAR `nwmlik.m`: `nwpost` on the dummies). Numbers (SumCoeff 0.1): varm -595.53 vs correct -630.38 (checked three ways: BEAR formula, general NW formula, and reproducing varm's value by switching the update off, diff 9e-11). Draws for FIXED hyper-parameters are unaffected (T4a' passes); the **hierarchical sampler for `soc`, `sur` (and lambda/alpha/psi jointly with them) uses this ML**, so their hyper-posteriors are suspect. (If this is deliberate following the GLP replication code, please state it; I could not run GLP's MATLAB code to confirm.)
4. Convention differences, not bugs: AR variance divisor (varm n-p vs BEAR n-2p-1), lag-decay exponent (alpha = 2*lambda3), unit-variance impact normalisation (Sigma/T vs Sigma/(T-k)), forecast bands (parameter uncertainty), strict vs weak sign checks, defaults.
5. BEAR-side observations (not varm bugs): `qzerores.m` order dependence; `nwprior_for_IV.m` names `invN0` what is actually the prior covariance and then inverts it (consistent but confusing); OLS-sign sampler draws Sigma with df T where T-k would be the diffuse-posterior df (varm mc uses T-m); BEAR comment `% T-k here!?`.

## Housekeeping / caveats

* All varm tests ran against the installed build in `rv/library/4.6/arm64` (built 00:08). Loading the working tree with `pkgload::load_all` triggered `compileAttributes()` and now fails to compile: `src/RcppExports.cpp` (modified 00:24) lacks `#include <RcppArmadillo.h>` that the HEAD version has (`error: use of undeclared identifier 'arma'`). Other sessions are editing the same tree (DESCRIPTION, R/*.R, NAMESPACE show as modified). I did not edit any varm R or C++ sources; my `load_all` may have been the trigger for the regeneration, so please check `git diff src/RcppExports.cpp` and re-run `Rcpp::compileAttributes()` in a session where RcppArmadillo is discoverable.
* Replica fidelity: `arloop`, `mprior`, `mpost`, `olsvar`, `olsforecast` replicas are verified against BEAR's stored numbers (B1-B6, A1-A5). `nwprior`, `nwpost`, `nwmlik`, `irfsim`, `fevd`, `hdecomp`, GK-IV replicas come from reading source only (nwprior/nwpost differ from mprior/mpost only by the Kronecker structure; T3 PASS and T11 PASS give mutual consistency with varm, which is an independent implementation).
* Monte-Carlo tolerances: |z| < 4.5 for means (several hundred comparisons per case); sd ratio within 3%.
* Reproduce: `cd varm && for f in validation/bear/0*.R; do Rscript $f; done` (needs python3 for the .mat reader, and `jsonlite`, `MASS`; `readxl` is taken from the site library).
