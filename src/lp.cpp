// [[Rcpp::depends(RcppArmadillo)]]
// Local projections (Jorda 2005): for each horizon h and response k
//   y_{k,t+h} = [S_t W_t]' theta + e,  theta = (beta (ns), gamma)
// OLS or 2SLS (instruments Z), Newey-West / EHW (HC1) standard errors for the first ns coefficients.
#include "varm_utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List cpp_lp(const arma::mat& Y, const arma::mat& S, const arma::mat& W, const arma::mat& Z, int H, int first,
            const arma::uvec& cum, int se_type, int nw_fixed, bool dof_adjust) {
  // se_type: 0 = Newey-West (Bartlett), 1 = HC1 (heteroskedasticity robust, no autocorrelation)
  const int T = Y.n_rows, K = Y.n_cols, ns = S.n_cols;
  const bool iv = Z.n_cols > 0;
  arma::mat cs = arma::cumsum(Y, 0);
  arma::cube beta(H + 1, K, ns), se(H + 1, K, ns);
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
    arma::mat G = arma::inv_sympd(Xh.t() * Xh);
    arma::mat Qm = Xh * G.cols(0, ns - 1);  // n x ns
    arma::mat Yh(n, K);
    for (int k = 0; k < K; ++k) {
      if (cum(k)) {
        for (int t = 0; t < n; ++t) {
          int tt = first + t;
          Yh(t, k) = cs(tt + h, k) - (tt > 0 ? cs(tt - 1, k) : 0.0);
        }
      } else Yh.col(k) = Y.rows(first + h, first + h + n - 1).col(k);
    }
    arma::mat Bh = G * (Xh.t() * Yh);       // m x K  (IV: (Xh'Xh)^-1 Xh' y)
    arma::mat E = Yh - X * Bh;
    int L = nw_fixed >= 0 ? nw_fixed : h + 1;
    double sc = dof_adjust || se_type == 1 ? (double)n / (n - m) : 1.0;
    for (int k = 0; k < K; ++k)
      for (int j = 0; j < ns; ++j) {
        beta(h, k, j) = Bh(j, k);
        arma::vec g = Qm.col(j) % E.col(k);
        double v = arma::dot(g, g);
        if (se_type == 0)
          for (int l = 1; l <= L && l < n; ++l) {
            double w = 1.0 - (double)l / (L + 1);
            v += 2 * w * arma::dot(g.subvec(l, n - 1), g.subvec(0, n - 1 - l));
          }
        se(h, k, j) = std::sqrt(sc * v);
      }
  }
  return List::create(_["beta"] = beta, _["se"] = se, _["nobs"] = nobs);
}
