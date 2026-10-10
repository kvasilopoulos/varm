# varm vs VAR-Toolbox (Cesa-Bianchi) - 2026-10-11

Reference: `external/VAR-Toolbox` (read-only; MATLAB/Octave not available and not installed). Method: read the `.m` source
(`VAR/VARmodel, VARmakexy, compute_wold/IR/VD/HD, recover_B, SR, SignRestrictions, VARdrawpost, LPmodel, OLSmodel`),
port the deterministic algorithms to R line by line (`validation/var-toolbox/mlref.R`), run both on the toolbox's bundled
data (xlsx converted to csv by `xlsx2csv.py`; BQ1989, Uhlig2005, GK2015, Primer) and compare.
Run: `Rscript validation/var-toolbox/run.R` from the package root (~1 min); log `run.log`, table `results.csv`.
No R sources edited, no commits. varm = working tree as of today (R 4.6.1; vars/svars/lpirfs not loadable, see baseline.md).

Caveat: the "MATLAB" side is my R port, validated only by reading the source, not by executing MATLAB. Numerical tolerances
are set from conditioning: the level-data VARs (Uhlig, GK) have cond(X'X) ~ 1e9-1e10, so normal-equation solutions agree
with each other to ~1e-6 only (see A).

Data/models: BQ = 2 vars (GDP growth, unemployment), p=8, T=159; Uhlig = 6 vars monthly, p=12, T=468; GK = 4 vars
(gs1, logcpi, logip, ebp) p=12, T=396, instrument FF4 observed in rows 127-396 (sentinel 123456789 -> NA); Primer = dlngdp, d(i1yr), T=123.
MATLAB layout IR(h,var,shock) is mapped to varm `[resp, shock, h]`; MATLAB step 1 = varm horizon 0.

## Summary
64 comparisons: 56 PASS, 8 FAIL. Every FAIL is a convention / documented-limitation difference, none is an arithmetic
disagreement. No outright varm bug found. Items worth a decision are marked (!) in "Discrepancies".

## Results
### A. Reduced-form estimation (`varm` vs `VARmodel`)
| test | how | number | status |
|---|---|---|---|
| coefficients, BQ p=8 | max abs diff, MATLAB column order [const, lag1..p] mapped to varm [lags, const] | 1.2e-14 | PASS |
| coefficients, Uhlig p=12 / GK p=12 | same; cond(X'X)=1e10 / 1e9 | 7.9e-7 / 5.3e-7 | PASS (tol 1e-5, see note) |
| vs QR solve (Uhlig/GK) | varm vs `qr.solve` | 3.0e-6 / 1.7e-6 | PASS (tol 1e-5) |
| ML normal-eq vs QR (Uhlig/GK) | the toolbox's own error | 3.8e-6 / 1.7e-6 | info |
| Sigma = RSS/(Te - ntotcoeff) | max diff, all 3 datasets | <= 2.3e-14 | PASS |
| residuals | max diff | 3e-14 (BQ), 1.7e-7 (Uhlig), 1.4e-8 (GK) | PASS |
| max companion eigenvalue | abs diff vs `stability()` | <= 1.1e-8 | PASS |
| const+trend (BQ p=4): lag coefs, trend slope, residuals | max diff | 3.5e-14, 1e-17, 2e-14 | PASS |
| const+trend: intercept | varm trend = 1..T, MATLAB trend = 1..Te (`VARmakexy`) | differs by exactly -p*slope (2.3e-3); residual after correcting 7e-16 | convention (FAIL raw, PASS adjusted) |

Note: varm and the toolbox both use normal equations (varm: `inv_sympd(X'X)`); on level data with p=12 both lose ~6 digits
vs QR. Harmless here; a QR/Cholesky-of-QR solve would be safer for long lags.

### B. Cholesky IRF (`irf(id_chol())` vs `compute_IR`, ident='short')
| test | number | status |
|---|---|---|
| GK 4 var, H=48, 1-sd, Wold recursion | 2.1e-8 (conditioning) | PASS |
| same vs companion-form recursion (`recurs='comp'`) | 2.1e-8 | PASS |
| unit shock (`impact=1`) vs `normalise=list(shock=,variable=,size=1)` | 1.2e-8 | PASS |
| `id_chol(order=c(3,1,4,2))` vs MATLAB on permuted data | 1.9e-9 | PASS |
| Uhlig 6 var, H=60 | 2.8e-7 | PASS |

### C. Long-run (`id_longrun()` vs ident='long'), BQ p=8, H=39
IRF all shocks 3.0e-14; impact matrix 3.0e-14; long-run matrix upper element 5e-16 (zero as required); BB'-Sigma 1.4e-17. **PASS**

### D. FEVD (`fevd` vs `compute_VD`, percent)
GK Cholesky H=48: 1.4e-5 pp (conditioning). BQ long-run H=39: 3.3e-12 pp. Shares sum to 1: 1e-16. **PASS**
(Same formula: cumulative squared IRF / cumulative MSE; MATLAB's MSE uses dof-corrected sigma, equal to sum of squared IRF over all shocks.)

### E. Historical decomposition (`hd` vs `compute_HD`)
Per-shock contributions, GK Cholesky: 1.7e-5; BQ long-run: 1.5e-13. varm `baseline` (actual - sum of shocks) vs MATLAB
HDinit+HDconst: 1.4e-5. **PASS.** Not tested: trend/exog/exoshock components (see NOT-TESTED).

### F. Proxy-SVAR / IV (`varm_proxy`, `id_proxy`, `proxy_B1` vs `recover_B` 'iv'), GK data
ML algorithm: drop NaN rows, first stage `OLSmodel(p,Z)` with constant, second stage on fitted p, sigma_b with
`1/(n-ntotcoeff)`, unit-variance scaling `sp`, sign from first-stage slope. varm: gamma = T^-1 sum u_t z_t (zero-filled), normalised to `size`.
| test | number | status |
|---|---|---|
| direction of gamma vs MATLAB Biv (cosine distance) | 1.2e-4 | PASS |
| normalised IRFs (gs1 impact = 1), H=48, zero-filled instrument vs ML NaN-drop | max 0.0211 abs = 1.6% of max IRF | FAIL (!) convention |
| same with instrument demeaned on its observed sample | 2.5e-8 | PASS |
| `id_proxy` vs `varm_proxy` point IRF | 5.6e-17 | PASS |
| unit-variance impact B1 vs MATLAB Biv | varm (0.2267, -5e-4, 5e-4, 0.1361) vs ML (0.1931, -4e-4, 5e-4, 0.1119): up to 17% | FAIL (!) convention |
| flipped-sign instrument leaves normalised IRF unchanged | 0 | PASS |
| FEVD of the proxy shock at h=48 (gs1) | ML 17.09% vs `fevd(id_proxy)` 100% | FAIL (!) documented limitation |
| FEVD recomputed from varm IRFs over total variance (4 vars) | ML 17.09/7.57/10.30/14.51; max diff 1.1 pp (from B1 scale) | PASS |
Why the IRF differs: ML centres through the first-stage constant on the instrument sample; varm's uncentred sum over the
zero-filled full sample equals the covariance only if the instrument is mean zero on its sample (FF4 mean -0.013, sd 0.053).
Why B1 differs: ML uses residual covariance on the instrument sample (n=270, dof n-ntotcoeff), varm uses the full-sample Sigma (T=384).
Bootstrap for IV (`mbb`, `wild`): code read, not executed against a port. varm MBB: block = round(5.03 T^.25) (ML floor), positionwise
recentring of residuals and proxy as in `VARmodel.m`; ML forbids iid `bs` with IV (varm has no `bs` for proxies, same stance). Block rule: 22 for both at T=384.

### G. Bootstrap bands (`draws(...,"bootstrap"/"wild")` vs `VARmodel` inference=1)
BQ long-run, 3000 draws each, 95% bands, percentile rule aside:
| test | number | status |
|---|---|---|
| mean band-width ratio varm/ML (BQ, iid bootstrap) | 1.075 | FAIL convention |
| GK Cholesky wild bootstrap, H=12, 600 draws | 1.063 | FAIL convention |
| ML port with varm's rescale `sqrt(Te/(Te-m))` = 1.0615 added: width ratio / draw-mean diff | 0.980 / 0.08 sd | PASS (explains it) |
| draw-mean IRF diff without rescale | 0.61 sd (cell h=0, y, shock 2) | FAIL, same cause |
| quantile rule, type 7 (varm) vs type 5 (MATLAB `prctile`), N=3000 | 0.029 sd | PASS |
| explosive draws (maxEig>=.9999) which ML rejects: kept by varm | 0.0% of 1000 (GK), 0.0% (BQ) | no effect on this data |

### H. Sign restrictions (Uhlig 2005 design: shock 1 = +ff, -pi, -comm, -nbres for 6 months; point estimate, 3000 accepted rotations)
| test | number | status |
|---|---|---|
| restrictions hold in 200 varm draws (horizons 0-5) | 0 violations | PASS |
| median IRF, max diff in units of rotation sd | 0.12 | PASS |
| 16/84% bands | 0.10 | PASS |
| 2.5/97.5% | 0.34 (tail MC noise) | PASS |
| acceptance rate | ML 0.351/rotation vs varm 0.0605 | info |
Both draw Haar rotations of chol(Sigma), so identified sets coincide. Acceptance differs by design: ML's 6x6 SIGN matrix
has 5 empty columns, so it scans all 6 rotated columns for shock 1 (per-column acceptance 0.058, any-column 0.35); varm
tests column 1 only. `accept_rate` is therefore not comparable across the two. ML check is `>=0` (zeros pass), varm `>0`.
Posterior (inference=1) path: ML draws (B, Sigma) via `VARdrawpost` per rotation search; varm `irf()` default is the bootstrap.

### I. Parameter-uncertainty draws (`VARdrawpost` vs `draws(,"mc")`, BQ, 4000 draws)
Sigma[1,1] mean: rel. diff 0.4%. Sd ratio varm/ML 1.045 (theory sqrt(T/(T-m)) = 1.062; MC noise ~2%). ML Wishart df = T,
varm df = T-m, so varm draws are more dispersed. **PASS** (convention).

### J. Local projections (`lp` vs `LPmodel`/`lp_single`)
Primer (dlngdp on d(i1yr), 4 lags, H=12, unit shock = ML `impact=1`): beta 3.7e-16 PASS. NW se: default varm bandwidth h+1 vs ML h:
1.7e-2 FAIL (convention); with `nw_lag = h`: 2.8e-16 PASS. LP-IV on GK (logip on gs1, instrument FF4, rows with FF4 only, 6 lags, H=12): beta 2.5e-10, se (`nw_lag=h`) 1.4e-9. **PASS**

### NOT-TESTED
- Narrative sign restrictions (`R.narr_sign`, `R.narr_dom`): no varm equivalent.
- `ident='exog'` (exogenous shock), `'sign+iv'` hybrid: no varm equivalent (partly emulated by `id_proxy` + sign, not compared).
- `VARopt.shut` (zeroing a variable's response), `longdiff` LP, ARDL, `nlag_ex` exogenous lags, `exoshock` HD components: no counterpart or not exercised.
- LP J-test (`jtest_chi2`), `impact=0` standard-deviation scaling in LP (varm is unit-shock; user rescales).
- Numerical IV bootstrap bands (mbb/wild) vs port; sign-restriction path with posterior draws end-to-end; HD for const=2.
- Plot functions (`VARirplot`, etc.).
- Anything needing MATLAB execution (rng streams are not comparable anyway; stochastic comparisons were distributional).

## Discrepancies
Suspected varm bugs: none confirmed. Items for the maintainer, ordered by importance:
1. (!) `fevd()` with `id_proxy()` returns 100% for the identified shock (denominator = identified column only). Documented in `id_proxy`,
   but the number is a misleading default; MATLAB reports 17.09% at h=48 (GK). Suggest total-variance denominators (Sigma-based MSE) for partial identification.
2. (!) `varm_proxy` / `id_proxy` doc says zero-filled missing instrument values leave normalised responses unchanged. Only exact for a mean-zero
   instrument on its sample; for FF4 (mean/sd = 0.25) the difference to the toolbox's NaN-drop+centred estimate is 1.6% of peak IRF. Either document, or
   demean on the observed sample (verified: then agreement is 2.5e-8). The unit-variance impact `B1` also differs by up to 17% because Sigma comes from the full sample, not the instrument sample.
3. Convention: bootstrap residuals are dof-rescaled by sqrt(Te/(Te-m)) and recentred; the toolbox uses raw OLS residuals. Bands are 6-7% wider
   (BQ p=8; more with long lags/short samples, e.g. 10% at m=49, Te=384). Not mentioned in `?draws`.
4. Convention: no rejection of explosive bootstrap draws (toolbox drops maxEig >= .9999). No effect on BQ/GK, may matter for near-unit-root data with long lags.
5. Convention: LP default Newey-West bandwidth h+1 (toolbox h); `nw_lag = h` reproduces the toolbox exactly.
6. Convention: trend regressor 1..T vs 1..Te (intercept differs by -p*slope, fit/residuals identical).
7. Convention: parameter draws use df T-m (toolbox T); percentile type 7 vs 5; sign-restriction acceptance counting; weak vs strict inequality.
8. Numerical: normal-equation OLS (also in the toolbox) loses ~6 digits on level data with p=12; harmless for IRFs (1e-5 level) but QR would remove it.
