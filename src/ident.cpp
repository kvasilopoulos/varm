// [[Rcpp::depends(RcppArmadillo)]]
#include "varm_utils.h"
using namespace Rcpp;

// ---------------------------------------------------------------- impact matrices
// Cholesky with optional variable ordering (ord is 0-based permutation)
static arma::mat impact_chol(const arma::mat& S, const arma::uvec& ord) {
  const int K = S.n_rows;
  arma::mat Sp = S.submat(ord, ord);
  arma::mat Lp = arma::chol(Sp, "lower");
  arma::mat P(K, K, arma::fill::zeros);
  P.submat(ord, ord) = Lp;
  return P;
}

// Blanchard-Quah: lower-triangular long-run multiplier. Ordering = variable order.
static arma::mat impact_lr(const arma::mat& B, const arma::mat& S, int K, int p) {
  arma::mat F = arma::eye(K, K);
  for (int l = 0; l < p; ++l) F -= B.rows(l * K, (l + 1) * K - 1).t();
  arma::mat C1 = arma::inv(F);
  arma::mat LR = arma::chol(C1 * S * C1.t(), "lower");
  return F * LR;
}

// Haar-distributed orthogonal matrix (QR of a Gaussian matrix, diag(R) > 0)
static arma::mat haar(int K) {
  arma::mat X(K, K);
  X.imbue([]() { return R::norm_rand(); });
  arma::mat Q, R;
  arma::qr(Q, R, X);
  for (int j = 0; j < K; ++j)
    if (R(j, j) < 0) Q.col(j) *= -1;
  return Q;
}

// Check sign restrictions S[i,j,h] in {-1,0,1} on shock column j. Returns +1 if column
// satisfies as is, -1 if its negative does, 0 if neither. Unrestricted column -> +1.
static int check_col(const arma::cube& Th, const arma::cube& S, int j) {
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

// ---------------------------------------------------------------- mapping
// Map reduced-form draws to structural IRFs.
// type: 0 chol, 1 long-run, 2 given P
// Returns cube K x K x (H+1)*N with slice index h + (H+1)*n  ->  R array [K,K,H+1,N]
// [[Rcpp::export]]
arma::cube cpp_map_irf(const arma::cube& Bs, const arma::cube& Ss, int p, int H, int type,
                       const arma::uvec& ord, const arma::mat& Pgiven) {
  const int N = Bs.n_slices, K = Ss.n_rows;
  arma::cube out(K, K, (H + 1) * N);
  for (int n = 0; n < N; ++n) {
    arma::mat P;
    if (type == 0) P = impact_chol(Ss.slice(n), ord);
    else if (type == 1) P = impact_lr(Bs.slice(n), Ss.slice(n), K, p);
    else P = Pgiven;
    arma::cube Th = varm_ma(Bs.slice(n), K, p, H);
    for (int h = 0; h <= H; ++h) out.slice(h + (H + 1) * n) = Th.slice(h) * P;
  }
  return out;
}

// Sign-restriction sampler (Rubio-Ramirez, Waggoner & Zha 2010; Uhlig 2005 uniform rotation).
// For each (B, Sigma) draw, try up to maxtry Haar rotations, keep up to nrot accepted.
// Returns accepted IRFs (K x K x (H+1)*Nacc), originating draw index (1-based) and impact matrices.
// [[Rcpp::export]]
List cpp_sign_irf(const arma::cube& Bs, const arma::cube& Ss, int p, int H, const arma::cube& S,
                  int nrot, int maxtry) {
  const int N = Bs.n_slices, K = Ss.n_rows;
  std::vector<arma::cube> acc;
  std::vector<int> src;
  std::vector<arma::mat> Ps;
  long tries = 0;
  const int Hs = S.n_slices;
  for (int n = 0; n < N; ++n) {
    arma::mat L = arma::chol(Ss.slice(n), "lower");
    arma::cube Phi = varm_ma(Bs.slice(n), K, p, std::max(H, Hs - 1));
    int got = 0;
    for (int t = 0; t < maxtry && got < nrot; ++t) {
      ++tries;
      arma::mat P = L * haar(K);
      arma::cube Th(K, K, Hs);
      for (int h = 0; h < Hs; ++h) Th.slice(h) = Phi.slice(h) * P;
      bool ok = true;
      for (int j = 0; j < K && ok; ++j) {
        int c = check_col(Th, S, j);
        if (c == 0) ok = false;
        else if (c < 0) P.col(j) *= -1;
        if (ok && c < 0) for (int h = 0; h < Hs; ++h) Th.slice(h).col(j) *= -1;
      }
      if (!ok) continue;
      arma::cube full(K, K, H + 1);
      for (int h = 0; h <= H; ++h) full.slice(h) = Phi.slice(h) * P;
      acc.push_back(full);
      src.push_back(n + 1);
      Ps.push_back(P);
      ++got;
    }
  }
  const int Na = acc.size();
  arma::cube out(K, K, (H + 1) * Na);
  arma::cube Pout(K, K, Na);
  for (int a = 0; a < Na; ++a) {
    for (int h = 0; h <= H; ++h) out.slice(h + (H + 1) * a) = acc[a].slice(h);
    Pout.slice(a) = Ps[a];
  }
  return List::create(_["irf"] = out, _["draw"] = wrap(src), _["P"] = Pout, _["tries"] = (double)tries,
                      _["accepted"] = Na);
}

// Rescale shock j so that the impact response of variable i equals `size` (0-based i, j)
// [[Rcpp::export]]
arma::cube cpp_normalise(arma::cube irf, int H1, int i, int j, double size) {
  const int N = irf.n_slices / H1;
  for (int n = 0; n < N; ++n) {
    double s = size / irf(i, j, H1 * n);
    for (int h = 0; h < H1; ++h) irf.slice(H1 * n + h).col(j) *= s;
  }
  return irf;
}

// Pointwise quantiles across draws: input K x K x (H1*N) -> K x K x (H1*nprob)
// output slice index h + H1*q
// [[Rcpp::export]]
arma::cube cpp_quantiles(const arma::cube& irf, int H1, const arma::vec& probs) {
  const int K = irf.n_rows, N = irf.n_slices / H1, Q = probs.n_elem;
  arma::cube out(K, K, H1 * Q);
  arma::vec v(N);
  for (int h = 0; h < H1; ++h)
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i) {
        for (int n = 0; n < N; ++n) v[n] = irf(i, j, h + H1 * n);
        arma::vec s = arma::sort(v);
        for (int q = 0; q < Q; ++q) {
          double pos = probs[q] * (N - 1);
          int lo = (int)std::floor(pos), hi = std::min(N - 1, lo + 1);
          out(i, j, h + H1 * q) = s[lo] + (pos - lo) * (s[hi] - s[lo]);
        }
      }
  return out;
}

// FEVD for many draws: input K x K x (H1*N) structural IRFs -> same shape, shares in [0,1]
// [[Rcpp::export]]
arma::cube cpp_fevd_many(const arma::cube& irf, int H1) {
  const int K = irf.n_rows, N = irf.n_slices / H1;
  arma::cube out(K, K, H1 * N);
  for (int n = 0; n < N; ++n) {
    arma::mat cum(K, K, arma::fill::zeros);
    for (int h = 0; h < H1; ++h) {
      cum += arma::square(irf.slice(H1 * n + h));
      arma::vec tot = arma::sum(cum, 1);
      out.slice(H1 * n + h) = cum.each_col() / tot;
    }
  }
  return out;
}
