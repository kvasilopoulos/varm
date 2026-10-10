// [[Rcpp::depends(RcppArmadillo)]]
// Local projections with the full cross-horizon covariance of the shock coefficients.
// Same regression as cpp_lp (OLS / 2SLS); beta_{h,j} is stacked as index h*ns + j, and
//   Cov(beta_{h,j}, beta_{h',j'}) = sum_{l} w_l (sum_t q_{t,h,j} e_{t,h} q_{t-l,h',j'} e_{t-l,h'} + transpose)
// with q = Xh (Xh'Xh)^-1 (first ns columns), Bartlett weights w_l = 1 - l/(L+1), L = max(h,h')+1 by default
// (so the diagonal reproduces cpp_lp), or a common fixed L (positive semi-definite).
#include "varm_utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List cpp_lp_cov(const arma::mat& Y, const arma::mat& S, const arma::mat& W, const arma::mat& Z, int H, int first,
                const arma::uvec& cum, int se_type, int nw_fixed, bool dof_adjust) {
  const int T = Y.n_rows, K = Y.n_cols, ns = S.n_cols;
  const bool iv = Z.n_cols > 0;
  const int Tn = T - first, M = (H + 1) * ns;
  arma::mat cs = arma::cumsum(Y, 0);
  arma::cube beta(H + 1, K, ns);
  arma::cube G(Tn, M, K, arma::fill::zeros);
  arma::vec nobs(H + 1);
  for (int h = 0; h <= H; ++h) {
    const int n = T - first - h;
    if (n <= (int)(ns + W.n_cols + 2)) stop("not enough observations for horizon %d", h);
    nobs(h) = n;
    arma::mat X = arma::join_rows(S.rows(first, first + n - 1), W.rows(first, first + n - 1));
    const int m = X.n_cols;
    arma::mat Xh = X;
    if (iv) {
      arma::mat Zf = arma::join_rows(Z.rows(first, first + n - 1), W.rows(first, first + n - 1));
      Xh = Zf * arma::solve(Zf, X);
    }
    arma::mat Gi = arma::inv_sympd(Xh.t() * Xh);
    arma::mat Qm = Xh * Gi.cols(0, ns - 1);
    arma::mat Yh(n, K);
    for (int k = 0; k < K; ++k) {
      if (cum(k)) {
        for (int t = 0; t < n; ++t) {
          int tt = first + t;
          Yh(t, k) = cs(tt + h, k) - (tt > 0 ? cs(tt - 1, k) : 0.0);
        }
      } else Yh.col(k) = Y.rows(first + h, first + h + n - 1).col(k);
    }
    arma::mat Bh = Gi * (Xh.t() * Yh);
    arma::mat E = Yh - X * Bh;
    double sq = std::sqrt(dof_adjust || se_type == 1 ? (double)n / (n - m) : 1.0);
    for (int k = 0; k < K; ++k)
      for (int j = 0; j < ns; ++j) {
        beta(h, k, j) = Bh(j, k);
        G.slice(k).col(h * ns + j).rows(0, n - 1) = sq * (Qm.col(j) % E.col(k));
      }
  }
  arma::cube V(M, M, K);
  arma::mat Lm(M, M);
  for (int a = 0; a < M; ++a)
    for (int b = 0; b < M; ++b) Lm(a, b) = nw_fixed >= 0 ? nw_fixed : std::max(a / ns, b / ns) + 1;
  const int Lmax = (int)Lm.max();
  for (int k = 0; k < K; ++k) {
    const arma::mat& Gk = G.slice(k);
    arma::mat tot = Gk.t() * Gk;
    if (se_type == 0)
      for (int l = 1; l <= Lmax && l < Tn; ++l) {
        arma::mat C = Gk.rows(l, Tn - 1).t() * Gk.rows(0, Tn - 1 - l);
        arma::mat Wl = arma::clamp(1.0 - l / (Lm + 1.0), 0.0, 1.0);
        Wl.elem(arma::find(Lm < l)).zeros();
        arma::mat WC = Wl % C;
        tot += WC + WC.t();
      }
    V.slice(k) = tot;
  }
  return List::create(_["beta"] = beta, _["vcov"] = V, _["nobs"] = nobs);
}
