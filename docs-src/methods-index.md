# varm methods index

One row per method: a one-line definition, the reference, the R entry point (linked
to its help page) and the compiled function (in `src/`) where there is one. The
"Chapter" column links the vignette that defines the method and explains how to read
it. Only methods that exist in the code are listed. Links point to the pkgdown site
(<https://kvasilopoulos.github.io/varm>).

Notation. $y_t$ is a $K$-vector, $A_l$ are lag matrices, $B=[A_1,\dots,A_p,c]'$,
$u_t$ the reduced-form innovation with covariance $\Sigma$, $u_t=P\varepsilon_t$ with
$PP'=\Sigma$ and $\varepsilon_t$ the structural shocks. $\Phi_h$ are the MA
coefficients (impulse responses to $u_t$), $\Theta_h=\Phi_hP$ the structural
IRF. Vignette names: [var-theory](https://kvasilopoulos.github.io/varm/articles/var-theory.html),
[identification](https://kvasilopoulos.github.io/varm/articles/identification.html),
[bayesian-var](https://kvasilopoulos.github.io/varm/articles/bayesian-var.html),
[tvp-var](https://kvasilopoulos.github.io/varm/articles/tvp-var.html),
[local-projections](https://kvasilopoulos.github.io/varm/articles/local-projections.html),
[lp-advanced](https://kvasilopoulos.github.io/varm/articles/lp-advanced.html),
[choosing-a-method](https://kvasilopoulos.github.io/varm/articles/choosing-a-method.html),
[validation](https://kvasilopoulos.github.io/varm/articles/validation.html).

## Reduced-form VAR

| Method | Definition | Reference | Where | Chapter |
|---|---|---|---|---|
| OLS / Gaussian ML | $y_t=c+\sum_{l=1}^pA_ly_{t-l}+u_t$; $\hat B=(X'X)^{-1}X'Y$, $\hat\Sigma=\hat U'\hat U/(T_e-m)$ ($T_e$ effective observations, $m$ regressors per equation; the ML version $\hat U'\hat U/T_e$ is also stored) | Lutkepohl (2005) | [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.html); `cpp_var_ols` (`core.cpp`) | var-theory |
| Companion stability | all eigenvalues of the companion matrix $F$ lie inside the unit circle | Hamilton (1994) | [`stability()`](https://kvasilopoulos.github.io/varm/reference/stability.html); `cpp_stability` | var-theory |
| Lag matrices | extract $A_1,\dots,A_p$ from $B$ | - | [`lag_matrices()`](https://kvasilopoulos.github.io/varm/reference/lag_matrices.html) | var-theory |
| Lag selection | $IC(p)=\log\det\hat\Sigma_{ML}(p)+C_T\,n_p/T$, $n_p=pK^2+Kd$ ($d$ deterministic terms), $C_T=2,\ 2\log\log T,\ \log T$ for AIC, HQ, SC; log FPE $=K\log\frac{T+Kp+d}{T-Kp-d}+\log\det\hat\Sigma_{ML}$, all on a common sample of $T=T_{\text{total}}-p_{\max}$ observations | Lutkepohl (2005) | [`var_select()`](https://kvasilopoulos.github.io/varm/reference/var_select.html); `cpp_var_ic` | var-theory |
| MA coefficients, IRF | $\Phi_0=I,\ \Phi_h=\sum_{l=1}^{\min(h,p)}\Phi_{h-l}A_l$; $\Theta_h=\Phi_hP$ | Lutkepohl (2005) | [`irf()`](https://kvasilopoulos.github.io/varm/reference/dynamics.html); `cpp_ma`, `cpp_map_irf` | var-theory |
| FEVD | $\omega_{ij}(h)=\sum_{s\le h}\Theta_{s,ij}^2\big/\sum_{s\le h}\sum_k\Theta_{s,ik}^2$ | Lutkepohl (2005) | [`fevd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.html); `cpp_fevd_many` | var-theory |
| Historical decomposition | $y_{t}=\text{deterministic}+\sum_{s=0}^{t-1}\Theta_s\varepsilon_{t-s}$, contribution of each shock | Burbidge and Harrison (1985) | [`hd()`](https://kvasilopoulos.github.io/varm/reference/dynamics.html); `cpp_hd` | var-theory |
| Forecasts | $\hat y_{T+h}=\hat B'x_{T+h}$ recursively; density by simulation | Lutkepohl (2005) | [`forecast()`](https://kvasilopoulos.github.io/varm/reference/dynamics.html); `cpp_path`, `cpp_forecast_draws` | var-theory |

## Inference: reduced-form draws

All return a `varm_draws` container of $(B^{(s)},\Sigma^{(s)})$.

| Method | Definition | Reference | Where | Chapter |
|---|---|---|---|---|
| Recursive-design residual bootstrap | resample $\hat u_t$ iid, regenerate $y^*_t$ recursively, re-estimate | Runkle (1987); Kilian (1998) | [`draws(method = "bootstrap")`](https://kvasilopoulos.github.io/varm/reference/draws.html); `cpp_var_boot` | var-theory |
| Wild bootstrap | $u^*_t=w_t\hat u_t$, $w_t$ Rademacher or Mammen | Goncalves and Kilian (2004); Mammen (1993) | `draws(method = "wild" / "wild_mammen")` | var-theory |
| Bias correction | subtract the bootstrap bias from $\hat B$, with stationarity adjustment (bootstrap-after-bootstrap) | Kilian (1998) | `draws(bias_correct = TRUE)` | var-theory |
| Asymptotic Monte Carlo | $\mathrm{vec}\,B\sim N(\mathrm{vec}\,\hat B,\hat\Sigma\otimes(X'X)^{-1})$, $\Sigma\sim IW$ | Lutkepohl (2005) | `draws(method = "mc")`; `cpp_var_mc` | var-theory |
| Percentile bands, normalisation, cumulation | quantiles of the draws of $\Theta_h$; rescale to unit impact; cumulate over $h$ | Kilian and Lutkepohl (2017) | `irf(level, normalise, cumulative)`; `cpp_quantiles`, `cpp_normalise` | identification, var-theory |

## Structural identification

| Method | Definition | Reference | Where | Chapter |
|---|---|---|---|---|
| Recursive (Cholesky) | $P=\mathrm{chol}(\Sigma)$ lower triangular, in the chosen ordering | Sims (1980) | [`id_chol()`](https://kvasilopoulos.github.io/varm/reference/identification.html) | identification |
| Long-run | $(I-\sum A_l)^{-1}P$ lower triangular | Blanchard and Quah (1989) | [`id_longrun()`](https://kvasilopoulos.github.io/varm/reference/identification.html) | identification |
| Sign restrictions | $P=LQ$ with $Q$ Haar-distributed orthogonal; keep draws with $\mathrm{sign}(\Theta_h)$ as restricted | Uhlig (2005); Rubio-Ramirez, Waggoner and Zha (2010) | [`id_sign()`](https://kvasilopoulos.github.io/varm/reference/identification.html); `cpp_sign_irf` | identification |
| Median-target rotation | the single accepted rotation closest to the pointwise median response | Fry and Pagan (2011) | [`irf_fry_pagan()`](https://kvasilopoulos.github.io/varm/reference/irf_fry_pagan.html) | identification |
| User impact matrix | $P$ supplied by the user | - | [`id_given()`](https://kvasilopoulos.github.io/varm/reference/identification.html) | identification |
| Zero and sign restrictions | null-space draw of $q_j$ with $e_i'\Phi_hLq_j=0$, then sign check (unweighted) | Arias, Rubio-Ramirez and Waggoner (2018) | [`id_zerosign()`](https://kvasilopoulos.github.io/varm/reference/id_zerosign.html); `cpp_zs_irf` | identification |
| Zero (A-B) restrictions | zeros in $P$ and in the long-run matrix; just-identified patterns solved by sequential projections | Blanchard and Quah (1989) | [`id_restrictions()`](https://kvasilopoulos.github.io/varm/reference/id_restrictions.html) | identification |
| Heteroskedasticity, two regimes | $\Sigma_1=PP'$, $\Sigma_2=P\Lambda P'$; $P$ from the joint eigen-decomposition | Rigobon (2003) | [`id_hetero()`](https://kvasilopoulos.github.io/varm/reference/id_hetero.html) | identification |

## Proxy-SVAR (external instruments)

| Method | Definition | Reference | Where | Chapter |
|---|---|---|---|---|
| Proxy identification | $\gamma=E[u_tz_t]\propto$ impact vector of the target shock, normalised so that the impact on `target` equals `size` | Mertens and Ravn (2013); Stock and Watson (2012) | [`id_proxy()`](https://kvasilopoulos.github.io/varm/reference/id_proxy.html); `cpp_proxy_irf` | identification |
| First-stage strength | (HAC) $F$ of the `target` residual on the proxy | Montiel Olea, Stock and Watson (2021) | [`proxy_fstat()`](https://kvasilopoulos.github.io/varm/reference/proxy_fstat.html) | identification |
| Proxy bootstrap | joint resampling of residuals and proxy (moving block or wild); one $\gamma$ per draw | Jentsch and Lunsford (2019) | [`proxy_draws()`](https://kvasilopoulos.github.io/varm/reference/proxy_draws.html), [`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.html); `cpp_proxy_boot` | identification |
| Weak-instrument robust sets | confidence set of the normalised response from inverting a delta-method Wald test: a bounded interval, the complement of an interval, or the whole line, depending on instrument strength | Montiel Olea, Stock and Watson (2021) | [`varm_proxy_robust()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_robust.html) | identification |
| Limited-information Bayesian proxy | $m_t=\Lambda'u_t+e_t$ with flat prior on $(\Lambda,\Sigma_e)$ on top of `bvarm()` draws | after Caldara and Herbst (2019) | [`varm_proxy_bayes()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_bayes.html) | identification |
| Bayesian proxy with sign and reliability restrictions | joint posterior of $(B,\Sigma,\text{rotation})$ under exogeneity, relevance and signs | Arias, Rubio-Ramirez and Waggoner (2021) | [`varm_proxy_arw()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_arw.html) (alias `id_proxy_bayes`); `cpp_proxy_rot` | identification |
| Prior-robust Bayesian proxy | bounds of the IRF over the class of priors for the rotation, given the reduced-form posterior | Giacomini, Kitagawa and Read (2022) | [`varm_proxy_gkr()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy_gkr.html) | identification |

## Bayesian VAR

| Method | Definition | Reference | Where | Chapter |
|---|---|---|---|---|
| Conjugate NIW posterior | $\beta\mid\Sigma\sim N(b,\Sigma\otimes\Omega)$, $\Sigma\sim IW(\Psi,M+2)$; direct draws | Koop and Korobilis (2010) | [`bvarm()`](https://kvasilopoulos.github.io/varm/reference/bvarm.html); `cpp_niw_fixed` (`bvar_niw.cpp`) | bayesian-var |
| Minnesota prior | prior variance of lag $l$ of variable $j$ in equation $i$: $\lambda^2\,\psi_i/(l^{\alpha}\psi_j)$ | Litterman (1986); Giannone, Lenza and Primiceri (2015) | [`prior_minnesota()`](https://kvasilopoulos.github.io/varm/reference/prior_minnesota.html) | bayesian-var |
| Sum-of-coefficients, single-unit-root dummies | dummy observations that shrink towards unit roots and cointegration | Doan, Litterman and Sims (1984); Sims (1993) | `prior_minnesota(soc, sur)` | bayesian-var |
| Closed-form marginal likelihood | $p(Y\mid\theta)$ of the NIW model with dummy observations | Giannone, Lenza and Primiceri (2015) | `cpp_niw_logpost`; `$logml` | bayesian-var |
| Hierarchical hyperparameters | $p(\theta\mid Y)\propto p(Y\mid\theta)p(\theta)$ explored by an optimiser and adaptive random-walk Metropolis | Giannone, Lenza and Primiceri (2015) | `bvarm()`; `cpp_niw_hier` | bayesian-var |
| Triangular sampler | $A_0y_t=\ldots+\Lambda_t^{1/2}\eta_t$ with $A_0$ unit lower triangular; equations sampled one at a time | Carriero, Clark and Marcellino (2019); Cogley and Sargent (2005) | [`bvarm(prior = prior_triangular())`](https://kvasilopoulos.github.io/varm/reference/prior_triangular.html); `cpp_bvar_tri` (`bvar_tri.cpp`) | bayesian-var |
| SSVS | spike-and-slab on each coefficient with inclusion indicator | George, Sun and Ni (2008) | `prior_triangular("ssvs")` | bayesian-var |
| Horseshoe | $\beta_j\sim N(0,\lambda_j^2\tau^2)$, half-Cauchy $\lambda_j,\tau$ | Carvalho, Polson and Scott (2010); Makalic and Schmidt (2016) | `prior_triangular("horseshoe")` | bayesian-var |
| Stochastic volatility | $\log\lambda_{it}=\log\lambda_{i,t-1}+\sigma_h\eta_{it}$; KSC 10-component mixture and forward-filtering backward-sampling | Kim, Shephard and Chib (1998); Omori et al. (2007); Carter and Kohn (1994) | `bvarm(sv = TRUE)`, [`volatility()`](https://kvasilopoulos.github.io/varm/reference/volatility.html); `cpp_sv_geweke`, `cpp_forecast_sv` | bayesian-var |
| MCMC diagnostics | ESS (Geyer), split R-hat, Geweke mean-equality $z$ | Vehtari et al. (2021); Geweke (1992) | [`mcmc_diagnostics()`](https://kvasilopoulos.github.io/varm/reference/mcmc_diagnostics.html) | bayesian-var |

## Time-varying parameter VAR

| Method | Definition | Reference | Where | Tidy / plot | Chapter |
|---|---|---|---|---|---|
| TVP-VAR-SV | $y_t=X_t\beta_t+A_t^{-1}\Lambda_t^{1/2}\varepsilon_t$ with $\beta_t,\ a_t,\ \log\lambda_t$ random walks; Gibbs sampler | Primiceri (2005); Cogley and Sargent (2005) | [`tvpvar()`](https://kvasilopoulos.github.io/varm/reference/tvpvar.html); `cpp_tvp` (`tvp.cpp`) | [`tidy()`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.html) (`date, term, estimate, conf.low, conf.high`), `glance()`, `autoplot(type = "coef")` | tvp-var |
| Training-sample prior | prior scales from an OLS fit on the first `tau` observations | Primiceri (2005) | [`tvp_prior()`](https://kvasilopoulos.github.io/varm/reference/tvp_prior.html) | - | tvp-var |
| Time-varying IRF, FEVD, forecasts | $\Theta_{h,t}$ computed at date $t$ from the stored paths; `irf()` rejects `id` and warns for out-of-sample dates | Del Negro and Primiceri (2015) | [`irf()`, `fevd()`, `forecast()` on `tvpvar`](https://kvasilopoulos.github.io/varm/reference/tvp_dynamics.html); `cpp_tvp_irf`, `cpp_tvp_forecast` | `autoplot(type = "irf")` | tvp-var |
| Response paths, volatilities | response of one variable to one shock at fixed horizon across dates; $\lambda_t^{1/2}$ paths | - | [`irf_path()`](https://kvasilopoulos.github.io/varm/reference/irf_path.html), [`tvp_volatility()`](https://kvasilopoulos.github.io/varm/reference/tvp_volatility.html); `cpp_tvp_irf_path`, `cpp_tvp_sd` | `irf_path()` is a tibble (`date, estimate, conf.low, conf.high`); `autoplot(type = "vol" / "shock_vol")` ([`autoplot`](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.html)) | tvp-var |

## Local projections

| Method | Definition | Reference | Where | Tidy / plot | Chapter |
|---|---|---|---|---|---|
| Local projections | $y_{t+h}=\mu_h+\beta_hs_t+\gamma_h'w_t+e_{t+h}$, one regression per $h$ | Jorda (2005) | [`lp()`](https://kvasilopoulos.github.io/varm/reference/lp.html); `cpp_lp` (`lp.cpp`) | [`tidy()`](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.html), `glance()`, [`autoplot()`](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.html) | local-projections |
| LP-IV | 2SLS of $y_{t+h}$ on the shock instrumented by $z_t$ | Stock and Watson (2018); Ramey (2016) | `lp(iv = )` | as `lp()` | local-projections |
| Newey-West and HC1 standard errors | Bartlett kernel with lag $h+1$; $\frac{n}{n-k}$-scaled sandwich | Newey and West (1987); White (1980) | `lp(se = )` | `std.error` column of `tidy()` | local-projections |
| Lag augmentation | add one uninterpreted lag and use HC standard errors | Montiel Olea and Plagborg-Moller (2021) | `lp(lag_augment = TRUE)` | as `lp()` | local-projections |
| Cumulative responses | $\sum_{j=0}^hy_{t+j}$ as the dependent variable | Jorda (2005) | `lp(cumulative = )` | as `lp()` | local-projections |
| Joint covariance, sup-t bands | $\mathrm{Cov}(\hat\beta_h,\hat\beta_{h'})$ from stacked HAC scores; $\hat\beta_h\pm c\,\hat\sigma_h$ with $c$ the $(1-\alpha)$ quantile of $\max_h\lvert Z_h\rvert$ | Montiel Olea and Plagborg-Moller (2019) | [`lp_joint()`](https://kvasilopoulos.github.io/varm/reference/lp_joint.html), [`lp_supt()`](https://kvasilopoulos.github.io/varm/reference/lp_supt.html); `cpp_lp_cov` | `tidy(lp_supt(...))`: `band` = `"sup-t"` / `"pointwise"`, `autoplot()` | lp-advanced |
| State-dependent LP | $y_{t+h}=F(z_{t-1})(\alpha^H_h+\beta^H_hs_t+\dots)+(1-F(z_{t-1}))(\alpha^L_h+\beta^L_hs_t+\dots)$, logistic $F$ | Auerbach and Gorodnichenko (2012, 2013) | [`lp_state()`](https://kvasilopoulos.github.io/varm/reference/lp_state.html) | `tidy()` with `state` column, `glance()`, `autoplot(type = "states" / "diff")` | lp-advanced |
| Smooth LP | $\beta_h=\sum_{k=1}^Jb_kB_k(h)$ in a B-spline basis, difference penalty on $b$ with weight chosen by cross-validation | Barnichon and Brownlees (2019) | [`lp_smooth()`](https://kvasilopoulos.github.io/varm/reference/lp_smooth.html) | `tidy(raw = TRUE)` adds the unsmoothed LP, `glance()`, `autoplot()` | lp-advanced |
| Panel LP | LP with unit/time fixed effects (alternating projections) and clustered standard errors | Jorda (2005) | [`lp_panel()`](https://kvasilopoulos.github.io/varm/reference/lp_panel.html) | `tidy()`, `glance()`, `autoplot()` | lp-advanced |
| Direct LP forecasts | $\hat y_{T+h}=\hat\mu_h+\hat\gamma_h'w_T$ from regressions of $y_{t+h}$ on $w_t$ | Jorda (2005) | [`lp_forecast()`](https://kvasilopoulos.github.io/varm/reference/lp_forecast.html) | `tidy()` (`horizon, variable, estimate, std.error, level, conf.low, conf.high`), `autoplot()` | lp-advanced |
| LP versus VAR | difference between LP and VAR responses with the same lags, in LP standard errors | Plagborg-Moller and Wolf (2021); Li, Plagborg-Moller and Wolf (2024) | [`compare_lp_var()`](https://kvasilopoulos.github.io/varm/reference/compare_lp_var.html); Monte Carlo in vignette | tibble (`horizon, response, lp, lp_se, var, diff, z`) | local-projections |

## Data and plotting

| Item | Where |
|---|---|
| Quarterly US inflation, unemployment, interest rate 1953Q1-2006Q3 | [`yraw`](https://kvasilopoulos.github.io/varm/reference/yraw.html) |
| Tidy output (tibbles, broom column names): VAR coefficients, IRF, FEVD, forecasts, historical decompositions, BVAR, TVP-VAR, LP | [`tidy()`](https://kvasilopoulos.github.io/varm/reference/tidy.varm.html), [`tidy()` on IRF/FEVD](https://kvasilopoulos.github.io/varm/reference/tidy.varm_irf.html) (`as.data.frame()` and `as_tibble()` are aliases), [forecast](https://kvasilopoulos.github.io/varm/reference/tidy.varm_forecast.html), [hd](https://kvasilopoulos.github.io/varm/reference/tidy.varm_hd.html), [BVAR](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.html), [TVP](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.html), [LP](https://kvasilopoulos.github.io/varm/reference/tidy.lpm.html) |
| ggplot2 figures (`autoplot()`; `plot()` prints the same object) | [IRF / FEVD](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.html), [forecast](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_forecast.html), [hd](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_hd.html), [draws](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_draws.html), [MCMC diagnostics](https://kvasilopoulos.github.io/varm/reference/autoplot.mcmc_diagnostics.html), [LP](https://kvasilopoulos.github.io/varm/reference/autoplot.lpm.html), [BVAR](https://kvasilopoulos.github.io/varm/reference/tidy.bvarm.html), [TVP](https://kvasilopoulos.github.io/varm/reference/tidy.tvpvar.html) |
| Theme and palette | [`theme_varm()`, `varm_palette()`](https://kvasilopoulos.github.io/varm/reference/theme_varm.html) |
| Shock normalisation | [`normalise_to()`](https://kvasilopoulos.github.io/varm/reference/normalise_to.html) |

## Validation of each method

For the evidence behind each row (reference implementation, agreement, coverage
simulations and known failures) see the
[validation report](https://kvasilopoulos.github.io/varm/articles/validation.html)
and `validation/`.
