// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "varm_utils.h"
using namespace Rcpp;

// Proxy-SVAR rotation sampler shared by the Arias-Rubio-Ramirez-Waggoner type posterior
// (varm_proxy_arw / id_proxy_bayes) and the Giacomini-Kitagawa-Read robust-Bayes bounds
// (varm_proxy_gkr).
//
// Per reduced-form draw phi = (B, Sigma, D): Sigma = Str Str' (lower Cholesky), the proxy equation
// m_t = D y_t + ... gives Cov(m_t, eps_t) = D Str Q. With r proxies for the LAST r shocks, the
// exogeneity restrictions are D Str q_j = 0, j = 1..K-r. Following GKR (2022, Algorithm 1 Steps 2.1-2.3)
// rotations are drawn column by column from the unit sphere of the orthogonal complement of
// {rows of D Str, previous columns}; the sign of column i is fixed by (Str^-1 e_i)'q_i >= 0
// (A0 has a non-negative diagonal) and the draw is kept if S * (Phi_h Str Q) >= 0 for all
// restricted (i, j, h). Up to `nq` accepted rotations per draw are used for the bounds
// [l, u] of the response of every variable at every horizon to shock `js` (GKR Algorithm 2, 3
// when tgt >= 0: responses to a shock with impact `size` on variable tgt); the first `keep`
// accepted rotations are returned in full.
// [[Rcpp::export]]
List cpp_proxy_rot(const arma::cube& Bs, const arma::cube& Ss, const arma::cube& Ds, int p, int H,
                   const arma::cube& S, int nq, int keep, int maxtry, int js, int tgt, double size) {
  const int N = Bs.n_slices, K = Ss.n_rows, r = Ds.n_rows, Hs = S.n_slices;
  const int Hm = std::max(H, Hs - 1);
  arma::cube lo(H + 1, K, N), up(H + 1, K, N);
  lo.fill(NA_REAL);
  up.fill(NA_REAL);
  arma::ivec nacc(N, arma::fill::zeros);
  std::vector<arma::cube> kept;
  std::vector<int> src;
  long tries = 0;
  for (int n = 0; n < N; ++n) {
    arma::mat Str = arma::chol(Ss.slice(n), "lower");
    arma::mat Si = arma::inv(arma::trimatl(Str));
    arma::cube Phi = varm_ma(Bs.slice(n), K, p, Hm);
    arma::mat Dt = Ds.slice(n) * Str;  // r x K
    arma::mat VF = arma::orth(Dt.t()); // K x r
    if ((int)VF.n_cols != r) continue;  // relevance fails
    arma::mat Q(K, K);
    arma::vec z(K);
    int got = 0;
    for (int t = 0; t < maxtry && got < nq; ++t) {
      ++tries;
      bool ok = true;
      for (int i = 0; i < K && ok; ++i) {
        z.imbue([]() { return R::norm_rand(); });
        if (i < K - r) z -= VF * (VF.t() * z);
        for (int j = 0; j < i; ++j) z -= Q.col(j) * arma::dot(Q.col(j), z);
        double nz = arma::norm(z);
        if (nz < 1e-10) { ok = false; break; }
        z /= nz;
        double s = arma::dot(Si.col(i), z);
        if (s == 0) s = R::unif_rand() < 0.5 ? -1.0 : 1.0;
        Q.col(i) = s < 0 ? -z : z;
      }
      if (!ok) continue;
      arma::mat P = Str * Q;
      for (int h = 0; h < Hs && ok; ++h) {
        arma::mat Th = Phi.slice(h) * P;
        for (int j = 0; j < K && ok; ++j)
          for (int i = 0; i < K; ++i) {
            double s = S(i, j, h);
            if (s != 0 && s * Th(i, j) < 0) { ok = false; break; }
          }
      }
      if (!ok) continue;
      ++got;
      if (js >= 0) {
        arma::vec g = P.col(js);
        double sc = tgt >= 0 ? size / g[tgt] : 1.0;
        for (int h = 0; h <= H; ++h) {
          arma::vec resp = Phi.slice(h) * g * sc;
          for (int i = 0; i < K; ++i) {
            if (got == 1 || resp[i] < lo(h, i, n)) lo(h, i, n) = resp[i];
            if (got == 1 || resp[i] > up(h, i, n)) up(h, i, n) = resp[i];
          }
        }
      }
      if (got <= keep) {
        arma::cube full(K, K, H + 1);
        for (int h = 0; h <= H; ++h) full.slice(h) = Phi.slice(h) * P;
        kept.push_back(full);
        src.push_back(n + 1);
      }
    }
    nacc[n] = got;
  }
  const int Nk = kept.size();
  arma::cube out(K, K, (H + 1) * Nk);
  for (int a = 0; a < Nk; ++a)
    for (int h = 0; h <= H; ++h) out.slice(h + (H + 1) * a) = kept[a].slice(h);
  return List::create(_["lower"] = lo, _["upper"] = up, _["nacc"] = nacc, _["irf"] = out,
                      _["draw"] = wrap(src), _["tries"] = (double)tries);
}
