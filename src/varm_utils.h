#ifndef VARM_UTILS_H
#define VARM_UTILS_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Coefficient layout used throughout varm:
//   y_t' = x_t' B + u_t',  x_t = (y_{t-1}', ..., y_{t-p}', d_t')'
//   B is m x K, m = K*p + d. A_l = B.rows((l-1)K, lK-1)'  (K x K).
//   Deterministics d_t are the LAST d rows (constant first, then trend).

// Regressor matrix for rows t = p..T-1 of Y (T x K): lags first, then det (T-p x d)
inline arma::mat varm_lag_design(const arma::mat& Y, int p, const arma::mat& det) {
  const arma::uword T = Y.n_rows, K = Y.n_cols, Te = T - p, d = det.n_cols;
  arma::mat X(Te, K * p + d);
  for (int l = 1; l <= p; ++l)
    X.cols((l - 1) * K, l * K - 1) = Y.rows(p - l, T - l - 1);
  if (d > 0) X.cols(K * p, K * p + d - 1) = det.rows(p, T - 1);
  return X;
}

// Companion matrix (Kp x Kp) from B (only the first Kp rows are used)
inline arma::mat varm_companion(const arma::mat& B, int K, int p) {
  arma::mat F(K * p, K * p, arma::fill::zeros);
  for (int l = 0; l < p; ++l) F.submat(0, l * K, K - 1, (l + 1) * K - 1) = B.rows(l * K, (l + 1) * K - 1).t();
  if (p > 1) F.submat(K, 0, K * p - 1, K * (p - 1) - 1) = arma::eye(K * (p - 1), K * (p - 1));
  return F;
}

// MA coefficients Phi_0..Phi_H (K x K x H+1) from B
inline arma::cube varm_ma(const arma::mat& B, int K, int p, int H) {
  arma::cube Phi(K, K, H + 1, arma::fill::zeros);
  Phi.slice(0) = arma::eye(K, K);
  for (int h = 1; h <= H; ++h)
    for (int l = 1; l <= std::min(h, p); ++l)
      Phi.slice(h) += B.rows((l - 1) * K, l * K - 1).t() * Phi.slice(h - l);
  return Phi;
}

// log multivariate gamma
inline double varm_lmgamma(double a, int n) {
  double r = n * (n - 1) / 4.0 * std::log(M_PI);
  for (int j = 0; j < n; ++j) r += R::lgammafn(a - j / 2.0);
  return r;
}

// Simulate y path from B given last p observations (oldest first), future deterministics and shocks
inline arma::mat varm_path(const arma::mat& B, int K, int p, const arma::mat& Ylast, const arma::mat& det,
                           const arma::mat& shocks) {
  const int H = shocks.n_rows, d = det.n_cols;
  arma::mat Y(p + H, K);
  Y.rows(0, p - 1) = Ylast;
  arma::rowvec x(K * p + d);
  for (int h = 0; h < H; ++h) {
    for (int l = 1; l <= p; ++l) x.subvec((l - 1) * K, l * K - 1) = Y.row(p + h - l);
    if (d > 0) x.subvec(K * p, K * p + d - 1) = det.row(h);
    Y.row(p + h) = x * B + shocks.row(h);
  }
  return Y.rows(p, p + H - 1);
}
#endif
