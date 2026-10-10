// [[Rcpp::depends(RcppArmadillo)]]
#include "varm_utils.h"
using namespace Rcpp;

// ---------------------------------------------------------------- OLS
// [[Rcpp::export]]
List cpp_var_ols(const arma::mat& Y, int p, const arma::mat& det) {
  const int K = Y.n_cols;
  arma::mat X = varm_lag_design(Y, p, det);
  arma::mat Ye = Y.rows(p, Y.n_rows - 1);
  const int Te = Ye.n_rows, m = X.n_cols;
  arma::mat XtXi;
  if (!arma::inv_sympd(XtXi, X.t() * X)) stop("X'X is singular: collinear regressors or too few observations");
  arma::mat B = XtXi * (X.t() * Ye);
  arma::mat U = Ye - X * B;
  arma::mat UU = U.t() * U;
  double ldet = arma::log_det_sympd(UU / Te);
  double ll = -0.5 * Te * K * std::log(2 * M_PI) - 0.5 * Te * ldet - 0.5 * Te * K;
  return List::create(_["B"] = B, _["U"] = U, _["Sigma"] = UU / (Te - m), _["Sigma_ml"] = UU / Te,
                      _["XtXinv"] = XtXi, _["X"] = X, _["Y"] = Ye, _["logLik"] = ll,
                      _["nobs"] = Te, _["npar"] = m);
}

// Information criteria for p = 1..pmax on a common sample (as in vars::VARselect)
// [[Rcpp::export]]
NumericMatrix cpp_var_ic(const arma::mat& Y, int pmax, const arma::mat& det) {
  const int K = Y.n_cols, d = det.n_cols;
  NumericMatrix out(pmax, 4);
  const int Te = Y.n_rows - pmax;
  arma::mat Ye = Y.rows(pmax, Y.n_rows - 1);
  for (int p = 1; p <= pmax; ++p) {
    arma::mat X(Te, K * p + d);
    for (int l = 1; l <= p; ++l) X.cols((l - 1) * K, l * K - 1) = Y.rows(pmax - l, Y.n_rows - l - 1);
    if (d > 0) X.cols(K * p, K * p + d - 1) = det.rows(pmax, Y.n_rows - 1);
    arma::mat U = Ye - X * arma::solve(X, Ye);
    double ld = arma::log_det_sympd(U.t() * U / Te);
    double np = p * K * K + K * d, T = Te;
    out(p - 1, 0) = ld + 2.0 / T * np;
    out(p - 1, 1) = ld + 2.0 * std::log(std::log(T)) / T * np;
    out(p - 1, 2) = ld + std::log(T) / T * np;
    out(p - 1, 3) = K * std::log((T + K * p + d) / (T - K * p - d)) + ld;  // log FPE (Lutkepohl)
  }
  return out;
}

// ---------------------------------------------------------------- dynamics
// [[Rcpp::export]]
arma::cube cpp_ma(const arma::mat& B, int K, int p, int H) { return varm_ma(B, K, p, H); }

// Structural IRFs Theta_h = Phi_h P, array [response, shock, h]
// [[Rcpp::export]]
arma::cube cpp_irf(const arma::mat& B, int K, int p, int H, const arma::mat& P) {
  arma::cube Phi = varm_ma(B, K, p, H);
  for (int h = 0; h <= H; ++h) Phi.slice(h) = Phi.slice(h) * P;
  return Phi;
}

// FEVD: share of the h-step forecast error variance of response i due to shock j
// (h index 0 = impact). Input: structural IRF cube Theta.
// [[Rcpp::export]]
arma::cube cpp_fevd(const arma::cube& Theta) {
  const int K = Theta.n_rows, H = Theta.n_slices;
  arma::cube out(K, K, H);
  arma::mat cum(K, K, arma::fill::zeros);
  for (int h = 0; h < H; ++h) {
    cum += arma::square(Theta.slice(h));
    arma::vec tot = arma::sum(cum, 1);
    out.slice(h) = cum.each_col() / tot;
  }
  return out;
}

// Eigenvalue moduli of the companion matrix
// [[Rcpp::export]]
arma::vec cpp_stability(const arma::mat& B, int K, int p) {
  arma::cx_vec ev = arma::eig_gen(varm_companion(B, K, p));
  return arma::sort(arma::abs(ev), "descend");
}

// Historical decomposition. Returns cube [t, i, j]: contribution of shock j
// to y_{t,i} (t = 1..Te) from the sample residuals; remainder = deterministic + initial conditions.
// [[Rcpp::export]]
arma::cube cpp_hd(const arma::mat& B, const arma::mat& U, const arma::mat& P, int K, int p) {
  const int Te = U.n_rows;
  arma::cube Theta = cpp_irf(B, K, p, Te - 1, P);
  arma::mat E = arma::solve(P, U.t()).t();  // structural shocks, Te x K
  arma::cube out(Te, K, K, arma::fill::zeros);
  for (int t = 0; t < Te; ++t)
    for (int s = 0; s <= t; ++s)
      for (int j = 0; j < K; ++j) {
        double e = E(t - s, j);
        for (int i = 0; i < K; ++i) out(t, i, j) += Theta(i, j, s) * e;
      }
  return out;
}

// ---------------------------------------------------------------- forecasting / simulation
// Point or simulated path. Ylast: p x K (oldest first, newest last).
// det: H x d future deterministic regressors. shocks: H x K innovations added each step (0 for point).
// [[Rcpp::export]]
arma::mat cpp_path(const arma::mat& B, int K, int p, const arma::mat& Ylast, const arma::mat& det,
                   const arma::mat& shocks) {
  return varm_path(B, K, p, Ylast, det, shocks);
}

// Density forecast over posterior/bootstrap draws: returns cube [draw, h, var]
// Bs: m x K x N; Ss: K x K x N. Parameter + shock uncertainty.
// [[Rcpp::export]]
arma::cube cpp_forecast_draws(const arma::cube& Bs, const arma::cube& Ss, int p, const arma::mat& Ylast,
                              const arma::mat& det, int H) {
  const int N = Bs.n_slices, K = Ylast.n_cols;
  arma::cube out(N, H, K);
  for (int n = 0; n < N; ++n) {
    arma::mat L = arma::chol(Ss.slice(n), "lower");
    arma::mat Z(H, K);
    Z.imbue([]() { return R::norm_rand(); });
    arma::mat path = varm_path(Bs.slice(n), K, p, Ylast, det, Z * L.t());
    for (int k = 0; k < K; ++k) out.slice(k).row(n) = path.col(k).t();
  }
  return out;
}

// ---------------------------------------------------------------- sampling of reduced-form uncertainty
// Residual / wild bootstrap of the reduced form. Returns draws of (B, Sigma).
// type: 0 = iid residual resampling, 1 = wild (Rademacher), 2 = wild (Mammen)
// bias: Kilian (1998) bootstrap-after-bootstrap correction of the point estimate
// Uses R's RNG (set.seed works).
// [[Rcpp::export]]
List cpp_var_boot(const arma::mat& Y, int p, const arma::mat& det, const arma::mat& B, const arma::mat& U,
                  int nboot, int type) {
  const int K = Y.n_cols, T = Y.n_rows, Te = T - p, d = det.n_cols;
  arma::cube Bs(B.n_rows, K, nboot), Ss(K, K, nboot);
  arma::mat Uc = U.each_row() - arma::mean(U, 0);
  Uc *= std::sqrt((double)Te / (Te - B.n_rows));  // dof rescale
  arma::mat Ystar(T, K);
  arma::rowvec x(K * p + d);
  const double s5 = std::sqrt(5.0), pm = (s5 + 1) / (2 * s5);
  for (int b = 0; b < nboot; ++b) {
    Ystar.rows(0, p - 1) = Y.rows(0, p - 1);
    arma::mat E(Te, K);
    if (type == 0) {
      for (int t = 0; t < Te; ++t) {
        int idx = std::min(Te - 1, (int)(R::unif_rand() * Te));
        E.row(t) = Uc.row(idx);
      }
    } else {
      for (int t = 0; t < Te; ++t) {
        double w;
        if (type == 1) w = R::unif_rand() < 0.5 ? -1.0 : 1.0;
        else w = R::unif_rand() < pm ? -(s5 - 1) / 2 : (s5 + 1) / 2;
        E.row(t) = Uc.row(t) * w;
      }
    }
    for (int t = p; t < T; ++t) {
      for (int l = 1; l <= p; ++l) x.subvec((l - 1) * K, l * K - 1) = Ystar.row(t - l);
      if (d > 0) x.subvec(K * p, K * p + d - 1) = det.row(t);
      Ystar.row(t) = x * B + E.row(t - p);
    }
    arma::mat X = varm_lag_design(Ystar, p, det);
    arma::mat Ye = Ystar.rows(p, T - 1);
    arma::mat Bb = arma::solve(X, Ye);
    arma::mat Ub = Ye - X * Bb;
    Bs.slice(b) = Bb;
    Ss.slice(b) = Ub.t() * Ub / (Te - B.n_rows);
  }
  return List::create(_["B"] = Bs, _["Sigma"] = Ss);
}

// Asymptotic Monte Carlo: B ~ N(Bhat, Sigma (x) (X'X)^-1), Sigma ~ IW(U'U, Te - m ...)
// [[Rcpp::export]]
List cpp_var_mc(const arma::mat& B, const arma::mat& U, const arma::mat& XtXinv, int ndraw) {
  const int m = B.n_rows, K = B.n_cols, Te = U.n_rows;
  arma::mat S = U.t() * U;
  arma::mat LX = arma::chol(XtXinv, "lower");
  arma::cube Bs(m, K, ndraw), Ss(K, K, ndraw);
  const double nu = Te - m;
  for (int n = 0; n < ndraw; ++n) {
    // Sigma ~ IW(S, nu): draw Wishart(S^-1, nu) via Bartlett
    arma::mat LSi = arma::chol(arma::inv_sympd(S), "lower");
    arma::mat A(K, K, arma::fill::zeros);
    for (int i = 0; i < K; ++i) {
      A(i, i) = std::sqrt(R::rchisq(nu - i));
      for (int j = 0; j < i; ++j) A(i, j) = R::norm_rand();
    }
    arma::mat LA = LSi * A;
    arma::mat W = LA * LA.t();
    arma::mat Sig = arma::inv_sympd(W);
    arma::mat LS = arma::chol(Sig, "lower");
    arma::mat Z(m, K);
    Z.imbue([]() { return R::norm_rand(); });
    Bs.slice(n) = B + LX * Z * LS.t();
    Ss.slice(n) = Sig;
  }
  return List::create(_["B"] = Bs, _["Sigma"] = Ss);
}
