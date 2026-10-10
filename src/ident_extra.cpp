// [[Rcpp::depends(RcppArmadillo)]]
#include "varm_utils.h"
using namespace Rcpp;

// Proxy-SVAR structural IRFs: shock column `sc` (0-based) = Theta_h g_n, other columns zero.
// G is K x Ng (Ng = 1 or N). Output K x K x (H+1)*N.
// [[Rcpp::export]]
arma::cube cpp_proxy_irf(const arma::cube& Bs, const arma::mat& G, int p, int H, int sc) {
  const int N = Bs.n_slices, K = G.n_rows;
  arma::cube out(K, K, (H + 1) * N, arma::fill::zeros);
  for (int n = 0; n < N; ++n) {
    arma::vec g = G.col(G.n_cols == 1 ? 0 : n);
    arma::cube Th = varm_ma(Bs.slice(n), K, p, H);
    for (int h = 0; h <= H; ++h) out.slice(h + (H + 1) * n).col(sc) = Th.slice(h) * g;
  }
  return out;
}

// Joint bootstrap of (residuals, proxy). type 0: wild (Rademacher weight shared by u_t and z_t,
// Jentsch & Lunsford 2019); type 1: moving block (blocks of length `block` of (u_t, z_t) pairs,
// Jentsch & Lunsford 2019 / Brüggemann et al. 2016). z has length Te. Returns B, Sigma, gamma (K x N).
// [[Rcpp::export]]
List cpp_proxy_boot(const arma::mat& Y, int p, const arma::mat& det, const arma::mat& B, const arma::mat& U,
                    const arma::vec& z, int nboot, int type, int block) {
  const int K = Y.n_cols, T = Y.n_rows, Te = T - p, d = det.n_cols, m = B.n_rows;
  arma::cube Bs(m, K, nboot), Ss(K, K, nboot);
  arma::mat Gs(K, nboot);
  arma::mat Uc = U.each_row() - arma::mean(U, 0);
  Uc *= std::sqrt((double)Te / (Te - m));
  arma::vec zc = type == 1 ? arma::vec(z - arma::mean(z)) : z;
  arma::mat Ystar(T, K), E(Te, K);
  arma::vec zb(Te);
  arma::rowvec x(K * p + d);
  block = std::max(1, std::min(block, Te));
  for (int b = 0; b < nboot; ++b) {
    if (type == 0) {
      for (int t = 0; t < Te; ++t) {
        double w = R::unif_rand() < 0.5 ? -1.0 : 1.0;
        E.row(t) = Uc.row(t) * w; zb[t] = zc[t] * w;
      }
    } else {
      int t = 0;
      while (t < Te) {
        int s = std::min(Te - block, (int)(R::unif_rand() * (Te - block + 1)));
        for (int k = 0; k < block && t < Te; ++k, ++t) { E.row(t) = Uc.row(s + k); zb[t] = zc[s + k]; }
      }
    }
    Ystar.rows(0, p - 1) = Y.rows(0, p - 1);
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
    Ss.slice(b) = Ub.t() * Ub / (Te - m);
    Gs.col(b) = Ub.t() * zb / Te;
  }
  return List::create(_["B"] = Bs, _["Sigma"] = Ss, _["gamma"] = Gs);
}

// ---------------------------------------------------------------------------------------
// Sign + zero restrictions (Arias, Rubio-Ramirez & Waggoner 2018, Algorithm 2 sampling steps).
static int check_col2(const arma::cube& Th, const arma::cube& S, int j) {
  bool pos = true, neg = true, any = false;
  for (arma::uword h = 0; h < S.n_slices; ++h)
    for (arma::uword i = 0; i < S.n_rows; ++i) {
      double s = S(i, j, h);
      if (s == 0) continue;
      any = true;
      double v = Th(i, j, h);
      if (s * v <= 0) pos = false;
      if (-s * v <= 0) neg = false;
      if (!pos && !neg) return 0;
    }
  if (!any) return 1;
  return pos ? 1 : -1;
}

// Z: K x K x Hz 0/1 array (1 = response of i to shock j at horizon h is zero); zlr: K x K 0/1 (long run).
// [[Rcpp::export]]
List cpp_zs_irf(const arma::cube& Bs, const arma::cube& Ss, int p, int H, const arma::cube& S,
                const arma::cube& Z, const arma::mat& zlr, int nrot, int maxtry) {
  const int N = Bs.n_slices, K = Ss.n_rows, Hs = S.n_slices, Hz = Z.n_slices;
  const int Hm = std::max(H, std::max(Hs - 1, Hz - 1));
  const bool haslr = arma::accu(zlr) > 0;
  // number of zero restrictions per shock, drawing order (most restricted first)
  arma::uvec nz(K, arma::fill::zeros);
  for (int j = 0; j < K; ++j) {
    for (int h = 0; h < Hz; ++h) for (int i = 0; i < K; ++i) if (Z(i, j, h) != 0) nz[j]++;
    for (int i = 0; i < K; ++i) if (zlr(i, j) != 0) nz[j]++;
  }
  arma::uvec ord = arma::stable_sort_index(nz, "descend");
  for (int k = 0; k < K; ++k)
    if ((int)nz[ord[k]] + k > K - 1)
      stop("zero restrictions violate the ARRW order condition: need #zeros(shock ordered k) <= K - k");
  std::vector<arma::cube> acc; std::vector<int> src; std::vector<arma::mat> Ps;
  long tries = 0;
  for (int n = 0; n < N; ++n) {
    arma::mat L = arma::chol(Ss.slice(n), "lower");
    arma::cube Phi = varm_ma(Bs.slice(n), K, p, Hm);
    arma::mat C1;
    if (haslr) {
      arma::mat F = arma::eye(K, K);
      for (int l = 0; l < p; ++l) F -= Bs.slice(n).rows(l * K, (l + 1) * K - 1).t();
      C1 = arma::inv(F);
    }
    std::vector<arma::mat> Zr(K);
    for (int j = 0; j < K; ++j) {
      arma::mat R(nz[j], K); int r = 0;
      for (int h = 0; h < Hz; ++h) for (int i = 0; i < K; ++i)
        if (Z(i, j, h) != 0) R.row(r++) = Phi.slice(h).row(i) * L;
      for (int i = 0; i < K; ++i) if (zlr(i, j) != 0) R.row(r++) = C1.row(i) * L;
      Zr[j] = R;
    }
    int got = 0;
    for (int t = 0; t < maxtry && got < nrot; ++t) {
      ++tries;
      arma::mat Q(K, K, arma::fill::zeros);
      bool ok = true;
      for (int k = 0; k < K && ok; ++k) {
        int j = ord[k];
        arma::mat V(nz[j] + k, K);
        if (nz[j] > 0) V.rows(0, nz[j] - 1) = Zr[j];
        for (int c = 0; c < k; ++c) V.row(nz[j] + c) = Q.col(ord[c]).t();
        arma::mat Nn;
        if (V.n_rows == 0) Nn = arma::eye(K, K); else Nn = arma::null(V);
        if (Nn.n_cols == 0) { ok = false; break; }
        arma::vec x(Nn.n_cols);
        x.imbue([]() { return R::norm_rand(); });
        Q.col(j) = Nn * x / arma::norm(x);
      }
      if (!ok) continue;
      arma::mat P = L * Q;
      arma::cube Th(K, K, Hs);
      for (int h = 0; h < Hs; ++h) Th.slice(h) = Phi.slice(h) * P;
      for (int j = 0; j < K && ok; ++j) {
        int c = check_col2(Th, S, j);
        if (c == 0) ok = false;
        else if (c < 0) { P.col(j) *= -1; for (int h = 0; h < Hs; ++h) Th.slice(h).col(j) *= -1; }
      }
      if (!ok) continue;
      arma::cube full(K, K, H + 1);
      for (int h = 0; h <= H; ++h) full.slice(h) = Phi.slice(h) * P;
      acc.push_back(full); src.push_back(n + 1); Ps.push_back(P); ++got;
    }
  }
  const int Na = acc.size();
  arma::cube out(K, K, (H + 1) * Na), Pout(K, K, Na);
  for (int a = 0; a < Na; ++a) {
    for (int h = 0; h <= H; ++h) out.slice(h + (H + 1) * a) = acc[a].slice(h);
    Pout.slice(a) = Ps[a];
  }
  return List::create(_["irf"] = out, _["draw"] = wrap(src), _["P"] = Pout, _["tries"] = (double)tries,
                      _["accepted"] = Na);
}
