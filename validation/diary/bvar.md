# Validation diary: varm vs naffe15/BVAR_ (MATLAB) — 2026-10-11

Reference: `external/BVAR_` (naffe15, `bvartools/`, v4.1, read-only). MATLAB/Octave are not installed, so every
overlapping deterministic algorithm was **ported line by line to R** (`validation/bvar/matlab_ports.R`, one function
per `.m` file: `rfvar3`, `varprior`, `matrictint`, `bvar_ml`, `iresponse`, `fevd`, `iresponse_longrun`,
`iresponse_proxy`, `forecasts`, `histdecomp`, `findQs/findP/generateDraw`, `ols_reg`, `IC`, plus a MATLAB v5 `.mat`
reader `readmat5.R`) and run on the toolbox's bundled data (`examples/BVAR tutorial/*.mat|xlsx`: `Data.mat`
euro-area monthly, `DataGK.mat` + `factor_data.xlsx` Gertler-Karadi, `DataGLP.mat`). Stochastic algorithms were compared
by distribution (two-sample z tests, KS tests, acceptance rates), deterministic ones to machine precision.
Caveat that applies to everything below: the ports were never run against real MATLAB output; they are checked only by
internal consistency (e.g. T1: BVAR_'s marginal likelihood equals the closed-form NIW marginal likelihood of its own
dummy-observation prior to 1e-9) and by exact identities.

Scripts (all in `validation/bvar/`, run with `sh run-all.sh`, output in `last-run.log`: **60 PASS, 17 FAIL**; FAIL lines are
classified below — most FAILs are defects/conventions of the *MATLAB* toolbox or convention differences, not varm bugs):
`01_minnesota.R 02_dynamics.R 03_proxy.R 04_signs.R 05_zerosign.R 06_hetero.R 07_lp.R 08_default_posterior.R`.

Environment: R 4.6.1; varm 0.1.0 at HEAD `d125cd9` + a dirty working tree (128 modified files; other agents were editing
the repo concurrently). I tested a **scratch copy** of the tree (`VARM_SRC`, loaded with `pkgload::load_all`) so that
compiling never touched the shared directory. Note: at 00:24 the tree's `src/RcppExports.cpp` had just been regenerated
without `#include <RcppArmadillo.h>` (compile error "undeclared identifier arma"); I re-ran
`Rcpp::compileAttributes()` in the scratch copy only. No file under `varm/R`, `varm/src`, nothing was edited, no git commit.

Status legend: PASS / FAIL (as printed by the scripts) / NOT-TESTED. "Attrib." says what a FAIL means:
**varm** = suspected varm issue, **conv** = convention difference (both defensible), **BVAR_** = defect/limitation of the MATLAB
toolbox, **num** = numerical-accuracy issue.

---------------------------------------------------------------------------------------------------------------------

## 1. Minnesota prior, marginal likelihood, posterior (`01_minnesota.R`)

### Mapping between the two parameterisations (derived and verified numerically)
BVAR_ uses Sims dummy observations with hyper-parameters (`tau, decay, lambda, mu, omega`); varm uses Giannone-Lenza-Primiceri
(GLP): `lambda, alpha, soc, sur, psi` with prior `Sigma ~ IW(diag(psi), M+2)`.

| BVAR_ | varm |
|---|---|
| prior var of lag-l coef of j in eq i = `Sigma_ii / (tau^2 l^(2 decay) sig_j^2)` | `lambda^2/(l^alpha psi_j)` => `alpha = 2 decay`, `lambda = sqrt(omega)/tau`, `psi = omega * sig^2` |
| `mu` (own persistence dummy `mu*diag(ybar)`) | `soc = 1/mu` |
| `lambda_B` (co-persistence, const column too) | `sur = 1/lambda_B` |
| `sig = std(y(1:lags,:))` (std of the first `lags` raw obs) | varm default is AR(p) residual variance; pass `psi` explicitly |
| `omega` = number of dummy rows on Sigma (df M*omega) | GLP fixes df = M+2. The two coincide only when `M*omega = M+3`, i.e. **M = 3, omega = 2** (BVAR_'s default omega = 2). All ML comparisons therefore use M = 3, omega = 2. |

I verified the mapping by comparing the full NIW prior implied by each code (B-prior, precision, IW scale, df): they agree
to 1e-7 (relative to the largest entry) once the mapping is applied.

### Tests
| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| T1 | R port of `bvar_ml.m` (rfvar3 SVD + `matrictint`) vs closed-form NIW marginal likelihood using BVAR_'s own prior parameters (Data.mat, 3 vars, p=2,4) | 3700.694283 vs 3700.694283; 3715.167565 vs 3715.167565 | PASS (port sanity) | — |
| T2b | varm `bvarm(prior_minnesota(..., hyper=character(0)))$logml` vs `bvar_ml` over 10 hyper grids (p=2,4; tau .3–10, decay .5–2, lambda 1–20, mu .5–5), BVAR_ ybar replaced by the rows varm uses (obs p+1..2p) | max abs diff 2.5e-5 (typically 1e-7) | PASS | — |
| T2a | same with BVAR_'s own `ybar = mean(y(1:p))` | diff −0.0071 … −0.0230, hyper-dependent | FAIL | **conv** (see D1) |
| T3 | varm logml insensitive to the constant's prior variance `var` (1e7 vs 1e12) when soc present | 1e-7 | PASS | — |
| T4 | M=2, omega=1: the Sigma priors differ (dof), so ML differs hyper-dependently | range −8.189 … −8.076 | PASS (documents the convention gap) | **conv** |
| T5a-c | posterior parameters at M=3, omega=2, ybar matched: B-hat, IW scale, IW dof | B 2.6e-8; S 3e-11 rel; nu 307 = 307 | PASS | — |
| T6 | 20 000 posterior draws: varm sampler vs BVAR_ sampler (`rand_inverse_wishart` + `chol(XXi) Z chol(Sigma)'`) — z-test of means of B and Sigma, variance ratio | max abs z 1.22 (B), 1.88 (Sigma); var ratio 0.971–1.015 | PASS | — |
| T7a | maximise ML over (tau, lambda) (decay=1, mu=2), ybar matched, Nelder-Mead on both surfaces | argmax (6.747, 2.307) = (6.747, 2.307); max logML 3707.944 = 3707.944 | PASS | — |
| T7b | same with BVAR_'s own ybar | argmax lambda 2.310 vs 2.307 (0.1%) | PASS | conv (tiny) |
| T8a/b | DataGLP (log GDP, C, I; p=5; BVAR_ defaults tau=3, decay=.5, lambda=5, mu=2, omega=2) | diff with BVAR_ ybar +0.0100; with varm's ybar rows 2.9e-8 | T8a FAIL / T8b PASS | **conv** (D1) |
| T9b | `soc=NULL` vs BVAR_ `mu=0`, constant made flat (`var=1e15`) | diffs ≤ 9e-5 | PASS | — |
| T9a | same with varm's default `var=1e7` | diffs 1.99, 4.88, 0.28, 7.92 log units | FAIL | **conv** (D4) |
| T10 | posterior means of Cholesky IRFs (4 horizons × 9 cells) from varm draws vs BVAR_ draws | max abs z 1.92 | PASS | — |

The hierarchical part (BVAR_ maximises the ML over log-hyperparameters with csminwel; varm maximises ML × gamma
hyper-priors with L-BFGS-B and then runs random-walk MH) is a different estimand; only the ML surface was compared
(T7). The hierarchical MH sampler itself has no BVAR_ counterpart: NOT-TESTED.

### Discrepancy notes
* **D1 (conv, low): `ybar` for the soc/sur dummies.** `src/bvar_niw.cpp::make_spec` takes `ybar = mean(Y.rows(0,p-1))` where
  `Y` is the *regression* sample (observations p+1..2p). BVAR_ (and Sims' rfvar3, whose doc says "first lags rows of ydata are
  initial conditions") uses the first p *raw* observations (1..p). It moves the log ML by 0.007–0.023 on Data.mat and 0.010 on
  DataGLP; the argmax of ML moves by 0.1%. varm is consistent with the CRAN `BVAR` package (its `bv_ml` also passes the trimmed
  `Y` to the dummy builders, and `tests/testthat/test-bvar-conj.R` matches it to 1e-6), so this is a convention, not an error;
  but using sample observations p+1..2p reuses likelihood information in the prior.
* **D4 (conv, medium for model comparison): `soc = NULL` makes the ML depend on `var`.** Without the soc rows only the sur row
  pins down the constant direction; the smallest eigenvalue of the prior precision is 3.5e-8 (BVAR_, flat constant) vs 1.3e-7
  (varm, var=1e7). With the default `var = 1e7`, `logml` differs from the flat-constant value by 0.3–7.9 log units on Data.mat
  (and 0.04 at var=1e9, 4e-4 at 1e11). Comparing ML across soc on/off, or against BVAR_, requires a much larger `var`.
* BVAR_ cannot evaluate `lambda = 0` (co-persistence off): the prior `XXi` is singular (my port raises in `chol`; MATLAB's
  `matrictint` only guards `diag(cx) < 100*eps`). varm handles `sur = NULL`.
* The prior on Sigma differs structurally (GLP `IW(Psi, M+2)` vs `M*omega` dummy rows with `Psi = omega*sig^2`), so BVAR_ with
  `omega` other than M=3/omega=2 cannot be reproduced by varm (T4).

---------------------------------------------------------------------------------------------------------------------

## 2. Reduced form, IRFs, FEVD, forecasts, historical decomposition, long-run (`02_dynamics.R`)
Data: DataGK `[logip logcpi gs1 ebp]`, T=396, lags=12, horizon 24 (example_3_irf.m case 1).

| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| P1 | varm OLS B vs QR-OLS (= BVAR_ rfvar3 without dummies) | 6.1e-8 (cond(X)=1.23e5) | PASS at 1e-6 | num |
| P1b | normal-equation residual max abs of X'(Y−XB) | varm 1.47, QR 8.9e-8 | FAIL | **num / varm** (D5) |
| P2 | `Sigma` = u'u/(T−k) (BVAR_ `Sigma_ols`) | 0 | PASS | — |
| D1 | Cholesky IRF, 1 std (`iresponse.m`) vs `irf(fit, id_chol(), ci="none")`, all 4×4×24 | 1.2e-15 | PASS | — |
| D2 | "1 percent" unit (`iresponse(...,unit=0)`) vs `normalise=list(shock=j,variable=j,size=1)` for each j | 4.9e-15 | PASS | — |
| D3 | FEVD at h=24 (`fevd.m`, percent) vs `fevd()` | 1.8e-14; rows sum to 100 | PASS | — |
| F1 | 12-step point forecast (`forecasts.m` no-shock path) vs `forecast(ci="none")` | 0 | PASS | — |
| F2a/b | forecast density at fixed (B,Sigma), 1e5 draws: mean = point path; covariance of h-step errors = Σ_i Θ_i Σ Θ_i' (h=1..6) | max|z| 2.36; max error 0.7% in correlation units (MC se 0.5%) | PASS | — |
| H1-H3 | historical decomposition (`histdecomp.m` at the OLS draw): contributions Te×K×K, BVAR_ adds-up check, deterministic+initial-condition part vs varm `baseline` | 3.8e-14, 1.3e-11, 1.3e-11 | PASS | — |
| H4 | HD with an orthogonal rotation Ω (`options.Omega`) vs `hd(fit, id_given(chol(Σ)Ω))` | 5.9e-14 | PASS | — |
| L1 | Blanchard-Quah (`iresponse_longrun.m`) vs `id_longrun()`, all 4 shocks × 24 h | 1.5e-15 | PASS | — |
| L3 | varm: long-run matrix C(1)P lower triangular | 1.2e-12 | PASS | — |
| L2a-c | BVAR_ `Q(1,1) = Q(1,1)*sign(Q(1,1))` flips only one *element*, not the column. Data.mat [IPI HICP Euribor1Y], p=1: raw Q(1,1) = −0.00744 | BVAR_ Q Q' − Σ = 1.4e-3 (not a factorisation); varm P P' − Σ = 3.5e-18; IRFs differ by 0.11 | L2a FAIL, L2b PASS, L2c FAIL | **BVAR_** (D6) |

* **D5 (num, low): varm solves OLS through `inv_sympd(X'X)`** (`src/core.cpp::cpp_var_ols`, also `cpp_lp`: `inv_sympd(Xh'Xh)`).
  With cond(X) = 1.2e5 (monthly levels, 12 lags) the coefficients are good to 6e-8 but the normal-equation residual is 1.5 (QR: 9e-8),
  and LP standard errors differ from a QR reference by ~1e-6 relative. Harmless for typical use; a QR/Cholesky-with-scaling
  would remove the cond² loss. All other comparisons in this section replace varm's `B, U, Sigma` by the QR values (engine-level
  comparison) so that conditioning does not mask logic: they then agree to 1e-14.
* **D6 (BVAR_ defect):** in 35 of 420 random variable orderings / lag orders over Data.mat, DataGK, DataGLP, DataEx the raw
  `Q(1,1)` is negative and BVAR_'s impact matrix stops satisfying Q Q' = Σ. varm has no such normalisation (sign of column 1 is
  then whatever the Cholesky of the long-run covariance gives; no corruption).

---------------------------------------------------------------------------------------------------------------------

## 3. Proxy-SVAR / external instrument (`03_proxy.R`)
Data: DataGK `[gs1 logip logcpi ebp]`, lags 12, hor 24; proxy `ff4_tc` (factor_data.xlsx col 4, 342 obs aligned to the last 342
residuals, example_3 case 6). BVAR_ `iresponse_proxy.m` = Mertens-Ravn/GK: 2SLS of the other residuals on the fitted first residual
(with intercepts, on the proxy sample), unit-variance scale `b11 = sqrt(Sigma11 − b12b12')`.

| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| X2 | `varm_proxy(fit, m = proxy centred on its sample, zero elsewhere, size = b11)` vs BVAR_ irs (24×4) | 6.7e-16 | PASS | — |
| X3 | `varm:::proxy_B1` (Mertens-Ravn B11/B21 formulae) with centred moments and df-corrected Sigma vs BVAR_ `b1` | 5.6e-17 | PASS | — |
| X5 | `irf(fit, id_proxy(z, fit))` vs `varm_proxy` (r=1) | 1.1e-16 | PASS | — |
| X6a | first-stage F (BVAR_ `F_m`) vs `proxy_fstat()` F (sub-sample with full proxy coverage) | 10.9652 = 10.9652 | PASS | — |
| X6b | robust F: BVAR_ White-HC0 Wald/n vs varm HC1 | 11.9977 vs 11.9275, ratio (T−2)/T exactly | PASS (small-sample factor) | conv |
| X7 | sub-sample with proxy covering every residual, raw (uncentred) proxy | 1.1e-12 | PASS | — |
| X1 | the actual GK case: varm zero-filled/uncentred proxy vs BVAR_ | IRF max abs diff 0.029; impact ratio 4.6% off; `varm_proxy$B1[1]` 0.2267 vs 0.2392 (5.2%) | FAIL | **conv / varm** (D7) |

* **D7 (conv, medium): zero-filling is not neutral.** `varm_proxy`/`id_proxy` use `gamma = T^-1 Σ u_t m_t` with `m` zero-filled
  outside the proxy sample. The docs claim zero-filling "leaves the normalised responses unchanged"; that is true only if the VAR
  residual mean over the proxy sample and the proxy mean are both zero. Over the proxy sample the residual mean is not zero
  (residuals have zero mean over the *full* sample only), and `ff4_tc` has a non-zero mean, so the effective statistic is
  Σ u m, not Cov(u, m). Centring `m` on its sample (X2) reproduces BVAR_/MR/GK exactly. Suggest centring `m` (and using the
  sample means of `u`) over the non-missing observations. Also `varm_proxy$B1` (unit-variance impact) uses the ML Sigma (÷T) whereas
  BVAR_ uses ÷(T−k): the unit-variance *scale* differs by about sqrt((T−k)/T) = 0.93 for T_e = 384, k = 49 (this, not centring, is
  most of the 5.2% gap in b11; the 4.6% gap in the impact *ratios* is the centring effect). With df-corrected Sigma and centred
  moments the two agree exactly (X3).

---------------------------------------------------------------------------------------------------------------------

## 4. Sign restrictions (`04_signs.R`)
At the OLS point estimate of DataGK `[logip logcpi gs1 ebp]`, p=12. Restrictions of example_3 case 2 (gs1 up and CPI down for
horizons 1–3 after shock 1). BVAR_: Haar rotation by `generateQ` (QR with sign fix), accept if all signs hold for `y` or
`−y` (whole array flipped, `checkrestrictions2.m`).

| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| S1a | acceptance probability, 60 000 BVAR_ rotations vs varm `id_sign` (6 000 accepted) | 0.3380 ± 0.0019 vs 0.3388 | PASS | — |
| S1b | accepted IRFs of shock 1, 4 vars × 12 horizons: two-sample z of means; quantile gaps (5,16,50,84,95%) | max abs z 1.41; max quantile gap 0.063 sd | PASS | — |
| S2 | restrictions on **two** shocks: BVAR_ flips the sign of the whole IRF array (one `fsign`), varm flips each column | acceptance 0.0882 (BVAR_) vs 0.1695 (varm), ratio 1.92 (theory 2) | PASS (ratio check) | **BVAR_** limitation (D8) |
| S3 | varm accepted impact matrices: P P' = Σ and sign restrictions | n=200, exact | PASS | — |

* **D8 (BVAR_ limitation):** with sign restrictions on more than one shock, BVAR_'s global sign flip accepts half as many
  (per shocks pair) rotations as the correct per-column treatment, and cannot realise e.g. (+,−) normalisations; varm's per-column
  flip enumerates the full sign-restricted set.

---------------------------------------------------------------------------------------------------------------------

## 5. Zero + sign restrictions (`05_zerosign.R`)
Impact restrictions of example_3 case 4 on DataGK (shock 3 has zeros on logip and logcpi; sign restrictions on shocks 1–3),
OLS point estimate. BVAR_ = `findQs`/`findP`/`generateDraw` (Binning 2013): random rotation `C = chol(Σ)Q0`, then a deterministic
orthogonal-complement construction; signs are *not* flipped, a column is rejected if its sign is wrong.

| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| Z1 | BVAR_ draws: zeros exact, W W' = Σ | 2e-16, 4e-16 | PASS | — |
| Z2 | varm `id_zerosign`: zeros exact, P P' = Σ (20 000 draws) | 0, 4e-16 | PASS | — |
| Z4a/b | exactly identified impact zeros (upper triangle): BVAR_ and `id_restrictions(impact=)` both give chol(Σ) up to column sign | 1e-16 / 0 | PASS | — |
| Z5a/b | exactly identified long-run zeros: both give the Blanchard-Quah matrix | 1e-15 / 4e-15 | PASS | — |
| Z6a/b | unfiltered law of the zero-restricted column q3: its angle in the 2-d null-space plane is uniform (KS) | varm D=0.0076 p=0.20; BVAR_ D=0.0066 p=0.35 | PASS / PASS | — |
| Z7a | varm: next column q1 \| q3 is conditional-Haar (coordinate in a fixed basis ~ U(−1,1)), KS | D=0.0059, p=0.49 | PASS | — |
| Z7b | BVAR_: same | D=0.048, p=7e-41 | FAIL | **BVAR_** (D9) |
| Z3, Z3b, Z3c | accepted impact matrices: varm vs BVAR_ (without / with column flipping): acceptance and 16-element distributions | acceptance 0.0161 (BVAR_, no flip), 0.1353 (BVAR_ flip) vs 0.1034 (varm); max abs z up to 10.7 / 23.3 | FAIL x3 | **BVAR_ / algorithm** (D9) |

* **D9:** both samplers satisfy the restrictions exactly and both draw the restricted column uniformly, but BVAR_'s remaining
  columns are deterministic functions of the first null-space vector and the random rotation, so the law over the remaining
  rotations is not conditionally Haar (Z7b); varm's sequential null-space sampler is (Z7a). Hence the sets/IRF bands differ
  (acceptance 0.135 vs 0.103). Neither carries the Arias-Rubio-Ramirez-Waggoner importance weights (varm documents this limitation
  in `?id_zerosign`; BVAR_ does not mention it). Not tested: long-run zero + sign combinations in the under-identified case.
  (Also by reading `findP.m`: `B(2:end,:)` assumes exactly one deterministic regressor — fails silently with a trend/exogenous block.)

---------------------------------------------------------------------------------------------------------------------

## 6. Identification through heteroskedasticity (`06_hetero.R`)
Simulated VAR(1), K=3, T=4000, known P0, variance ratios (4, 1.6, 0.5), regime = halves. BVAR_ `iresponse_heterosked.m`:
`X = S0\S1; [~,S,V] = svd(X); iresponse(Phi, S, hor, V)`.

| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| HT1 | varm `id_hetero` recovers P0 up to order/sign | max rel error 0.042 | PASS | — |
| HT2 | varm lambda estimates | 3.959, 1.707, 0.506 vs 4, 1.6, 0.5 | PASS | — |
| HT3 | varm formula on exact population covariances | 2.3e-16 | PASS | — |
| HT4 | BVAR_ impact `chol(diag(svals)) V` reproduces the regime-1 covariance | rel error 1.86 (W W' = diag(svals) by construction) | FAIL | **BVAR_** (D10) |
| HT5 | BVAR_ impact recovers P0 | no column match | FAIL | **BVAR_** |
| HT6 | BVAR_ formula on exact population covariances | not equal to P0 | FAIL | **BVAR_** |

* **D10 (BVAR_ suspected bug, by reading + port):** the code passes the *singular value matrix* `S` of `inv(S0)*S1` as the
  covariance to `iresponse` and rotates by the right singular vectors `V`, so the impact matrix is `diag(sqrt(s)) V`, which satisfies
  W W' = diag(s) rather than Σ and is not the Rigobon matrix (correct: `chol(S0)` times eigenvectors of the symmetrised
  `L0^-1 S1 L0^-T`, which is what varm does). Not run in MATLAB; the conclusion rests on the algebra and the port.
  varm also uses uncentred crossproducts per regime (BVAR_: `cov`, centred) — immaterial for VAR residuals with a constant.

---------------------------------------------------------------------------------------------------------------------

## 7. Local projections (`07_lp.R`)
DataGK `[logip logcpi gs1 ebp]`, lags 4, horizon 12. BVAR_ `directmethods.m` runs OLS LPs of `y_{t+h}` on `[lags(y), controls, proxy, const]`,
standard errors Hamilton/Newey-West (`robust_se_=1`, `L = lags+h+1`); `ir_lp` is the Cholesky response `beta_{h-1}(lag-1 block)' chol(Sigma_0)`.

| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| L1a | varm `lp(shock=gs1, contemp=logip,logcpi)` beta vs QR-OLS `ols_reg` on the same design, h=0..12, all responses | max abs diff / se 2.5e-8 | PASS | — |
| L1b/c | NW se (varm `nw_lag=h+1`, and `nw_lag=lags+h+1` = BVAR_'s bandwidth) vs `ols_reg.m` Hamilton formula | relative 1.6e-6, 2.6e-6 | PASS | num (D5) |
| L2a | BVAR_ LP-Cholesky responses (`beta_{h-1}` lag-1 coefficients × Ω) vs varm unit-impact LP with contemporaneous controls, h≥1 | ≤ 0.81 se_varm per shock (0.81, 0.24, 0.43, 0.21) | PASS (<1 se) | conv — different regressions (not an algebraic identity); also differs from the VAR-Cholesky IRF by up to 1.8 |
| L3b | varm proxy LP engine (`lp(shock=numeric proxy)`) equals QR-OLS/NW on its own design | 1.4e-6 rel | PASS | — |
| L3 | BVAR_ proxy LP (`[lags, proxy, const]`) vs varm (adds `lags` lags of the external shock), h ∈ {0,3,6,12} | max 0.71 se | FAIL (threshold 0.5 se) | **conv**: varm always controls for lags of an external shock; BVAR_ does not |
| — | confidence bands | BVAR_ uses t(T−ny−1)=1.6488, varm z=1.6449 | note | conv |

BVAR_'s Bayesian LP (`blp_ml`, conjugate prior with `tau` shrinkage, `p2p.m`) has no varm counterpart: NOT-TESTED.

---------------------------------------------------------------------------------------------------------------------

## 8. Default posterior and information criteria (`08_default_posterior.R`)
| ID | What / how | Numbers | Status | Attrib. |
|---|---|---|---|---|
| D1 | BVAR_ default sampler: E[Sigma] = S/(df−K−1), df = T_e − k − (K+1) | max abs z 1.26 | PASS | — |
| D2 | varm `draws(fit,"mc")`: E[Sigma] = U'U/(nu−K−1), nu = T_e − k | max abs z 1.20 | PASS | — |
| D3 | do the two Sigma posteriors coincide? | BVAR_ df 380, varm nu 385; E[Sigma] ratio 0.98684 | FAIL | **conv** (D11) |
| D4 | varm mc: Cov(vec B) = E[Sigma] ⊗ (X'X)^-1 | max rel dev 2.8% (N=4e4; 0.6% at N=2e5) | PASS | — |
| I1 | lag choice from BVAR_ `IC.m` vs varm `var_select`, 1..6 lags (DataGK) | BVAR_ 6,6,6 vs varm 6,3,2 (AIC,HQ,BIC) | FAIL | **BVAR_/conv** (D12) |
| I2 | varm `var_select` vs `vars::VARselect` | 6,3,2 = 6,3,2 | PASS | — |

* **D11 (conv):** BVAR_ with no Minnesota prior sets `flat = 1`, i.e. posterior df = T−k−(K+1): IW(S, T−k−K−1), the posterior
  of a flat prior on Σ, although the code/docs call it Jeffreys (Jeffreys gives IW(S, T−k), which is what `cpp_var_mc` samples).
  Difference: E[Sigma] 1.3% on this example.
* **D12:** `IC.m` evaluates the log-likelihood with the df-corrected Σ and penalises with `K` = regressors per equation (not
  K·N + N(N+1)/2 parameters) and a different sample for each lag, so BIC/HQ keep selecting the maximum lag. varm matches
  `vars::VARselect`.

---------------------------------------------------------------------------------------------------------------------

## 9. Summary of findings

### Suspected varm issues (all minor/moderate; no wrong IRF/FEVD/HD/forecast/LP/ML formula found)
1. **(num)** `cpp_var_ols` / `cpp_lp` use normal equations; with ill-conditioned designs (cond 1e5) coefficients are good to ~6e-8 and
   X'e residuals are far from orthogonal (1.5 vs 1e-7). Prefer QR. (D5)
2. **(conv→possible bias)** `varm_proxy`/`id_proxy` zero-fill missing proxy values and use uncentred moments: on the GK data
   the impact ratios differ 4.6% (IRF max 0.029) from the centred Mertens-Ravn/GK estimator; the doc claim that zero-filling "leaves the
   normalised responses unchanged" does not hold in general. Centre `m` and `u` on the observed proxy sample. (D7)
3. **(conv)** soc/sur dummies use `ybar` from observations p+1..2p instead of the first p raw observations (shifts ML by ~0.01);
   same as CRAN `BVAR`. (D1)
4. **(conv/usability)** with `soc = NULL`, `logml` depends strongly on `var` (default 1e7 is 0.3–8 log units from the flat-constant
   limit); document or default larger. (D4)
5. LP with an external shock always adds lags of the shock to the controls (unlike BVAR_); LP bands use normal quantiles. (L3)

### Defects found in BVAR_ (the reference), for the record
* `iresponse_longrun.m`: `Q(1,1) = Q(1,1)*sign(Q(1,1))` flips a single element; wrong impact matrix whenever Q(1,1)<0 (~8% of cases tried). (D6)
* `iresponse_heterosked.m`: uses singular values as covariance and V as rotation; not an identification (W W' = diag(s)). (D10)
* `checkrestrictions2.m`/`iresponse_sign.m`: single global sign flip; wrong/overly restrictive with multi-shock restrictions. (D8)
* `findP.m` sampler: non-conditional-Haar draws of the unrestricted columns; assumes exactly one deterministic term. (D9)
* No-Minnesota posterior labelled Jeffreys but is IW(S, T−k−K−1). `IC.m` penalty too small. `lambda = 0` breaks the ML. (D11, D12)

### Features where varm and BVAR_ agree to machine precision / statistically
reduced-form OLS and Sigma; Cholesky IRF (1 std and unit-impact); FEVD; point forecasts and forecast-error covariance; historical
decomposition (also rotated); Blanchard-Quah (when BVAR_'s normalisation is not triggered); exact zero-restriction identification
(impact and long-run); sign-restriction acceptance rate and accepted-IRF distribution (single shock); proxy-SVAR IRFs/F-stat given
matching conventions; Minnesota NIW marginal likelihood (exactly, after mapping hyper-parameters, ybar rows and Sigma-prior df); NIW posterior B-hat/S/dof and draw moments;
LP point estimates and Newey-West standard errors on identical designs.

## 10. NOT-TESTED (and why)
* In BVAR_ but **absent from varm**: narrative sign restrictions (`iresponse_sign_narrative`), higher-moment/ICA identification
  (`iresponse_sign_hmoments*`), Pandemic Priors, Priors for the Long Run (PLR), covid/heteroskedastic weights, FAVAR, mixed-frequency
  VAR/Kalman filter, `bdfm_`, nowcasting, connectedness, Ridge/Lasso/Elastic-net VAR, pooled panel VAR, VARX with exogenous block,
  robust-to-kurtosis/skewness Bayes, conditional forecasts (`cforecasts*`), Bayesian LP (`blp_ml`), Giacomini-Kitagawa robust credible regions
  via rotations (varm has `varm_proxy_gkr`, a different construction), the trend-cycle filters (`cmintools`, `Hpfilter`, …).
* In varm but **absent from BVAR_**: TVP-VAR, stochastic volatility, SSVS/horseshoe/triangular priors, bootstrap/wild-bootstrap
  inference, AR/MOSW robust proxy intervals, ARW/GKR Bayesian proxy SVAR, `lp_state/smooth/joint/supt/panel`, Fry-Pagan.
* Hierarchical hyper-parameter MH sampler vs BVAR_'s `bvar_max_hyper` (different estimands; only the ML surface compared, T7).
* Sign restrictions with posterior (Bayesian) reduced-form draws end-to-end and FEVD under set identification: only point-estimate
  and acceptance-set equivalence were tested.
* Anything requiring actual MATLAB execution (e.g. exact random-number streams, `csminwel` optimiser trajectories, plotting).
