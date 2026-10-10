// [[Rcpp::depends(RcppArmadillo)]]
// TVP-VAR with stochastic volatility (Primiceri 2005; Del Negro & Primiceri 2015 ordering):
//   y_t = X_t' B_t + A_t^{-1} Lambda_t^{1/2} eps_t,   beta_t = vec(B_t) = beta_{t-1} + eta_t,  eta ~ N(0,Q)
//   alpha_t (free elements of unit-lower-triangular A_t) = alpha_{t-1} + zeta_t, zeta_i ~ N(0,S_i) (blocks)
//   h_it = log lambda_it = h_{i,t-1} + nu_it, nu ~ N(0,W_i)             (KSC 10-component mixture)
// Gibbs blocks, in order: beta | A,h,Q ; alpha | beta,h,S ; (s,h) | beta,alpha ; Q, S, W.
// State draws use the Durbin-Koopman (2002) simulation smoother: O(T * ns^2 * ny) instead of
// O(T * ns^3) for Carter-Kohn, and no Cholesky of near-singular conditional state covariances.
#include "varm_utils.h"
using namespace Rcpp;

namespace {

const double MIX_P[10] = {0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115};
const double MIX_M[10] = {1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000};
const double MIX_V[10] = {0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342};

inline double rig(double a, double b) { return 1.0 / R::rgamma(a, 1.0 / b); }

inline arma::vec rnorm_vec(int n) {
  arma::vec z(n);
  for (int i = 0; i < n; ++i) z(i) = R::norm_rand();
  return z;
}

// Lower Cholesky factor, with eigenvalue clipping fallback (for numerically semi-definite input)
arma::mat chol_safe(const arma::mat& S) {
  arma::mat L;
  if (arma::chol(L, S, "lower")) return L;
  arma::vec ev; arma::mat V;
  arma::eig_sym(ev, V, 0.5 * (S + S.t()));
  ev = arma::sqrt(arma::clamp(ev, 0.0, arma::datum::inf));
  return V * arma::diagmat(ev);
}

// Inverse-Wishart(Psi, nu) draw via Bartlett decomposition of the Wishart(nu, Psi^{-1})
arma::mat riwish(const arma::mat& Psi, double nu) {
  const int n = Psi.n_rows;
  arma::mat L = chol_safe(arma::inv_sympd(0.5 * (Psi + Psi.t())));
  arma::mat Ab(n, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    Ab(i, i) = std::sqrt(R::rchisq(nu - i));
    for (int j = 0; j < i; ++j) Ab(i, j) = R::norm_rand();
  }
  arma::mat LA = L * Ab;
  return arma::inv_sympd(LA * LA.t());
}

// KSC log-volatility update (see bvar_tri.cpp). e: orthogonalised residuals, h: path, W: innovation variance.
void sv_update(const arma::vec& e, arma::vec& h, double& h0, double& W, double mu0, double V0, double a_w,
               double b_w, std::vector<int>& s) {
  const int T = e.n_elem;
  arma::vec ys(T);
  for (int t = 0; t < T; ++t) ys(t) = std::log(e(t) * e(t) + 1e-4);
  double lp[10];
  for (int t = 0; t < T; ++t) {
    double mx = -1e300;
    for (int k = 0; k < 10; ++k) {
      double r = ys(t) - h(t) - MIX_M[k];
      lp[k] = std::log(MIX_P[k]) - 0.5 * std::log(MIX_V[k]) - 0.5 * r * r / MIX_V[k];
      mx = std::max(mx, lp[k]);
    }
    double tot = 0;
    for (int k = 0; k < 10; ++k) { lp[k] = std::exp(lp[k] - mx); tot += lp[k]; }
    double u = R::unif_rand() * tot, c = 0;
    int kk = 9;
    for (int k = 0; k < 10; ++k) { c += lp[k]; if (u <= c) { kk = k; break; } }
    s[t] = kk;
  }
  arma::vec a(T + 1), P(T + 1);
  a(0) = mu0; P(0) = V0;
  for (int t = 0; t < T; ++t) {
    double Pm = P(t) + W, v = MIX_V[s[t]], F = Pm + v, Kg = Pm / F;
    a(t + 1) = a(t) + Kg * (ys(t) - MIX_M[s[t]] - a(t));
    P(t + 1) = Pm * (1.0 - Kg);
  }
  h(T - 1) = a(T) + std::sqrt(P(T)) * R::norm_rand();
  for (int t = T - 2; t >= -1; --t) {
    double At = a(t + 1), Pt = P(t + 1);
    double nxt = (t == T - 2) ? h(T - 1) : h(t + 1);
    double g = Pt / (Pt + W);
    double mean = At + g * (nxt - At), var = Pt * W / (Pt + W);
    double draw = mean + std::sqrt(var) * R::norm_rand();
    if (t >= 0) h(t) = draw; else h0 = draw;
  }
  double ss = (h(0) - h0) * (h(0) - h0);
  for (int t = 1; t < T; ++t) ss += (h(t) - h(t - 1)) * (h(t) - h(t - 1));
  W = rig(a_w + 0.5 * T, b_w + 0.5 * ss);
}

// Durbin-Koopman simulation smoother for
//   y_t = Z_t s_t + e_t, e_t ~ N(0, R_t);   s_1 ~ N(b0, P0), s_t = s_{t-1} + eta_t, eta ~ N(0, Q)
// Y: ny x T, Z: ny x ns x T, R: ny x ny x T. P0c, Qc: lower Cholesky factors. Returns ns x T draw.
arma::mat dk_draw(const arma::cube& Z, const arma::cube& R, const arma::mat& Y, const arma::vec& b0,
                  const arma::mat& P0, const arma::mat& P0c, const arma::mat& Q, const arma::mat& Qc) {
  const int ny = Y.n_rows, T = Y.n_cols, ns = b0.n_elem;
  // 1. simulate (s+, y+) and form y* = y - Z b0 - y+
  arma::mat sp(ns, T), ystar(ny, T);
  arma::vec g = P0c * rnorm_vec(ns);
  for (int t = 0; t < T; ++t) {
    if (t > 0) g += Qc * rnorm_vec(ns);
    sp.col(t) = g;
    arma::mat Rc = chol_safe(R.slice(t));
    arma::vec yp = Z.slice(t) * g + Rc * rnorm_vec(ny);
    ystar.col(t) = Y.col(t) - Z.slice(t) * b0 - yp;
  }
  // 2. Kalman filter on y* (zero prior mean)
  arma::vec a(ns, arma::fill::zeros);
  arma::mat P = P0;
  arma::cube Kg(ns, ny, T);
  arma::mat Fiv(ny, T);
  for (int t = 0; t < T; ++t) {
    const arma::mat& Zt = Z.slice(t);
    arma::mat PZ = P * Zt.t();
    arma::mat F = Zt * PZ + R.slice(t);
    arma::mat Fi = arma::inv_sympd(0.5 * (F + F.t()));
    arma::vec v = ystar.col(t) - Zt * a;
    Kg.slice(t) = PZ * Fi;
    Fiv.col(t) = Fi * v;
    a += Kg.slice(t) * v;
    P -= Kg.slice(t) * PZ.t();
    P = 0.5 * (P + P.t()) + Q;
  }
  // 3. backward recursion for r_t, then smoothed state mean
  arma::vec r(ns, arma::fill::zeros);
  arma::mat rr(ns, T);
  for (int t = T - 1; t >= 0; --t) {
    rr.col(t) = r;
    const arma::mat& Zt = Z.slice(t);
    r = Zt.t() * (Fiv.col(t) - Kg.slice(t).t() * r) + r;
  }
  arma::mat out(ns, T);
  arma::vec sh = P0 * r;                      // smoothed state at t = 1
  for (int t = 0; t < T; ++t) {
    if (t > 0) sh += Q * rr.col(t - 1);
    out.col(t) = b0 + sp.col(t) + sh;
  }
  return out;
}

// build A^{-1} (K x K) from alpha column (unit lower triangular, row-wise packing)
inline arma::mat build_Ainv(const arma::vec& al, int K) {
  arma::mat A(K, K, arma::fill::eye);
  int k = 0;
  for (int i = 1; i < K; ++i)
    for (int j = 0; j < i; ++j) A(i, j) = al(k++);
  return arma::inv(arma::trimatl(A));
}

inline arma::mat build_A(const arma::vec& al, int K) {
  arma::mat A(K, K, arma::fill::eye);
  int k = 0;
  for (int i = 1; i < K; ++i)
    for (int j = 0; j < i; ++j) A(i, j) = al(k++);
  return A;
}

// largest root of the VAR(p) companion for coefficient vector beta (m x K layout)
inline bool explosive(const arma::vec& beta, int m, int K, int p) {
  arma::mat B = arma::reshape(beta, m, K);
  arma::cx_vec ev = arma::eig_gen(varm_companion(B, K, p));
  return arma::abs(ev).max() >= 1.0;
}

}  // namespace

// Gibbs sampler. Y, X: estimation sample (T x K, T x m). Prior pieces are computed in R from the training sample.
// opts: [sv, tv_beta, tv_A, q_diag, reject, max_tries]
// [[Rcpp::export]]
List cpp_tvp(const arma::mat& Y, const arma::mat& X, int p, const arma::vec& b0, const arma::mat& P0,
             const arma::vec& a0, const List& Pa0, const arma::mat& Qbar, double nuQ, const arma::mat& Q0,
             const List& Sbar, const arma::vec& nuS, const List& S0, const arma::vec& hmu, double hV0,
             const arma::vec& wprior, const arma::vec& s2prior, const arma::vec& opts, int ndraw, int nburn,
             int thin, int nsp) {
  const int T = Y.n_rows, K = Y.n_cols, m = X.n_cols, nb = m * K, nal = K * (K - 1) / 2;
  const bool sv = opts(0) > 0.5, tvb = opts(1) > 0.5, tva = opts(2) > 0.5, qdiag = opts(3) > 0.5,
             reject = opts(4) > 0.5;
  const int max_tries = (int)opts(5);

  // observation matrices Z_t = I_K (x) x_t'  (K x nb)
  arma::cube Zb(K, nb, T, arma::fill::zeros);
  for (int t = 0; t < T; ++t)
    for (int i = 0; i < K; ++i) Zb.slice(t).submat(i, i * m, i, (i + 1) * m - 1) = X.row(t);
  arma::mat Yt = Y.t();
  arma::mat P0inv = arma::inv_sympd(P0), P0c = chol_safe(P0);
  std::vector<arma::mat> Pa0inv(K), Pa0c(K), Pa0m(K), S(K), Sbi(K);
  std::vector<int> aidx(K);
  for (int i = 1; i < K; ++i) {
    arma::mat Pi = as<arma::mat>(Pa0[i]); Pa0m[i] = Pi;
    Pa0inv[i] = arma::inv_sympd(Pi); Pa0c[i] = chol_safe(Pi);
    S[i] = as<arma::mat>(S0[i]); Sbi[i] = as<arma::mat>(Sbar[i]);
    aidx[i] = i * (i - 1) / 2;
  }

  // state
  arma::mat Bt(nb, T), At(std::max(nal, 1), T, arma::fill::zeros), H(K, T);
  Bt.each_col() = b0;
  for (int t = 0; t < T; ++t) At.col(t).head(nal) = a0;
  arma::vec h0 = hmu, W(K), s2(K);
  for (int i = 0; i < K; ++i) { H.row(i).fill(hmu(i)); W(i) = wprior(1) / wprior(0); s2(i) = std::exp(hmu(i)); }
  arma::mat Q = Q0;
  std::vector<int> smix(T);

  // storage
  arma::cube Bp(nb, T, nsp), Ap(std::max(nal, 1), T, nsp), Hp(K, T, nsp), Qp(nb, nb, nsp);
  arma::mat Sp(std::max(1, (int)(K * (K - 1) * (2 * K - 1) / 6)), nsp, arma::fill::zeros), Wp(K, nsp);
  arma::mat bsum(nb, T, arma::fill::zeros), asum(std::max(nal, 1), T, arma::fill::zeros), hsum(K, T, arma::fill::zeros);
  arma::mat qsd(ndraw, nb), ssd(ndraw, std::max(nal, 1)), wsd(ndraw, K);
  const int spstep = std::max(1, (ndraw + nsp - 1) / std::max(nsp, 1));
  int saved = 0, stored = 0;
  const int total = nburn + ndraw * thin;
  arma::mat U(K, T);
  arma::cube Rc(K, K, T);

  for (int it = 0; it < total; ++it) {
    // Ainv_t and Sigma_t = Ainv Lambda Ainv'
    std::vector<arma::mat> Ainvs(T);
    for (int t = 0; t < T; ++t) {
      Ainvs[t] = build_Ainv(At.col(t), K);
      arma::mat AL = Ainvs[t].each_row() % arma::exp(0.5 * H.col(t)).t();
      Rc.slice(t) = AL * AL.t();
    }
    // 1. beta
    if (tvb) {
      arma::mat Qc = chol_safe(Q);
      arma::mat Bnew;
      int tries = 0;
      while (true) {
        Bnew = dk_draw(Zb, Rc, Yt, b0, P0, P0c, Q, Qc);
        ++tries;
        bool bad = false;
        if (reject) for (int t = 0; t < T && !bad; ++t) bad = explosive(Bnew.col(t), m, K, p);
        if (!bad) { Bt = Bnew; break; }
        if (tries >= max_tries) break;  // keep previous path
      }
    } else {
      arma::mat prec = P0inv;
      arma::vec rhs = P0inv * b0;
      for (int t = 0; t < T; ++t) {
        arma::mat Si = build_A(At.col(t), K);        // A_t
        arma::mat Sinv = Si.t() * arma::diagmat(arma::exp(-H.col(t))) * Si;
        arma::mat ZS = Zb.slice(t).t() * Sinv;
        prec += ZS * Zb.slice(t);
        rhs += ZS * Yt.col(t);
      }
      arma::mat L = chol_safe(arma::inv_sympd(0.5 * (prec + prec.t())));
      arma::vec mu = arma::solve(prec, rhs);
      arma::vec bn = mu + L * rnorm_vec(nb);
      Bt.each_col() = bn;
    }
    // residuals U = y - X'B_t
    for (int t = 0; t < T; ++t) U.col(t) = Yt.col(t) - Zb.slice(t) * Bt.col(t);
    // 2. alpha (equation by equation)
    for (int i = 1; i < K; ++i) {
      const int ni = i;
      arma::cube Zc(1, ni, T);
      arma::cube Rcc(1, 1, T);
      arma::mat Yi(1, T);
      for (int t = 0; t < T; ++t) {
        Zc.slice(t) = -U.submat(0, t, ni - 1, t).t();
        Rcc(0, 0, t) = std::exp(H(i, t));
        Yi(0, t) = U(i, t);
      }
      arma::vec a0i = a0.subvec(aidx[i], aidx[i] + ni - 1);
      if (tva) {
        arma::mat Sci = chol_safe(S[i]);
        arma::mat path = dk_draw(Zc, Rcc, Yi, a0i, Pa0m[i], Pa0c[i], S[i], Sci);
        for (int t = 0; t < T; ++t) At.col(t).subvec(aidx[i], aidx[i] + ni - 1) = path.col(t);
      } else {
        arma::mat prec = Pa0inv[i];
        arma::vec rhs = Pa0inv[i] * a0i;
        for (int t = 0; t < T; ++t) {
          prec += Zc.slice(t).t() * Zc.slice(t) / Rcc(0, 0, t);
          rhs += Zc.slice(t).t() * Yi.col(t) / Rcc(0, 0, t);
        }
        arma::mat L = chol_safe(arma::inv_sympd(0.5 * (prec + prec.t())));
        arma::vec an = arma::solve(prec, rhs) + L * rnorm_vec(ni);
        for (int t = 0; t < T; ++t) At.col(t).subvec(aidx[i], aidx[i] + ni - 1) = an;
      }
    }
    // 3. volatilities
    for (int i = 0; i < K; ++i) {
      arma::vec e(T);
      for (int t = 0; t < T; ++t) {
        double v = U(i, t);
        for (int j = 0; j < i; ++j) v += At(aidx[i] + j, t) * U(j, t);
        e(t) = v;
      }
      if (sv) {
        arma::vec h = H.row(i).t();
        sv_update(e, h, h0(i), W(i), hmu(i), hV0, wprior(0), wprior(1), smix);
        H.row(i) = h.t();
      } else {
        s2(i) = rig(s2prior(0) + 0.5 * T, s2prior(1) * std::exp(hmu(i)) + 0.5 * arma::dot(e, e));
        H.row(i).fill(std::log(s2(i)));
      }
    }
    // 4. innovation variances
    if (tvb) {
      arma::mat D = arma::diff(Bt, 1, 1);          // nb x (T-1)
      arma::mat SS = D * D.t();
      if (qdiag) {
        for (int j = 0; j < nb; ++j)
          Q(j, j) = rig(0.5 * nuQ + 0.5 * (T - 1), 0.5 * Qbar(j, j) + 0.5 * SS(j, j));
      } else {
        Q = riwish(Qbar + SS, nuQ + T - 1);
      }
    }
    if (tva) {
      for (int i = 1; i < K; ++i) {
        arma::mat D = arma::diff(At.rows(aidx[i], aidx[i] + i - 1), 1, 1);
        S[i] = riwish(Sbi[i] + D * D.t(), nuS(i) + T - 1);
      }
    }
    // storage
    if (it >= nburn && ((it - nburn) % thin == thin - 1)) {
      bsum += Bt; asum += At; hsum += H;
      for (int j = 0; j < nb; ++j) qsd(saved, j) = tvb ? std::sqrt(Q(j, j)) : 0.0;
      for (int i = 1; i < K; ++i)
        for (int j = 0; j < i; ++j) ssd(saved, aidx[i] + j) = tva ? std::sqrt(S[i](j, j)) : 0.0;
      for (int i = 0; i < K; ++i) wsd(saved, i) = sv ? std::sqrt(W(i)) : 0.0;
      if (((saved + 1) % spstep == 0) && stored < nsp) {
        Bp.slice(stored) = Bt; Ap.slice(stored) = At; Hp.slice(stored) = H; Qp.slice(stored) = tvb ? Q : arma::mat(nb, nb, arma::fill::zeros);
        int off = 0;
        for (int i = 1; i < K; ++i) {
          arma::vec sv_ = tva ? arma::vec(arma::vectorise(S[i])) : arma::vec(i * i, arma::fill::zeros);
          Sp.col(stored).subvec(off, off + i * i - 1) = sv_;
          off += i * i;
        }
        for (int i = 0; i < K; ++i) Wp(i, stored) = sv ? W(i) : 0.0;
        ++stored;
      }
      ++saved;
    }
    if ((it & 15) == 0) Rcpp::checkUserInterrupt();
  }
  const int nst = stored;
  return List::create(_["beta"] = arma::cube(Bp.slices(0, nst - 1)), _["alpha"] = arma::cube(Ap.slices(0, nst - 1)),
                      _["h"] = arma::cube(Hp.slices(0, nst - 1)), _["Q"] = arma::cube(Qp.slices(0, nst - 1)),
                      _["S"] = arma::mat(Sp.cols(0, nst - 1)), _["W"] = arma::mat(Wp.cols(0, nst - 1)),
                      _["beta_mean"] = bsum / ndraw, _["alpha_mean"] = asum / ndraw, _["h_mean"] = hsum / ndraw,
                      _["q_sd"] = qsd, _["s_sd"] = ssd, _["w_sd"] = wsd);
}

// Frozen-parameter IRFs at selected dates. Returns K x K x ((H+1) * N * nd): slice = ((d * N + n) * (H+1) + h),
// element (i, j) = response of variable i to a one-s.d. structural shock j: Phi_h(B_t) A_t^{-1} Lambda_t^{1/2}.
// [[Rcpp::export]]
arma::cube cpp_tvp_irf(const arma::cube& Bp, const arma::cube& Ap, const arma::cube& Hp, int m, int K, int p, int H,
                       const arma::uvec& dates) {
  const int N = Bp.n_slices, nd = dates.n_elem;
  arma::cube out(K, K, (H + 1) * N * nd);
  for (int d = 0; d < nd; ++d)
    for (int n = 0; n < N; ++n) {
      const int t = dates(d);
      arma::mat B(const_cast<double*>(Bp.slice(n).colptr(t)), m, K, false);
      arma::cube Phi = varm_ma(B, K, p, H);
      arma::mat Ainv = build_Ainv(Ap.slice(n).col(t), K);
      arma::mat P = Ainv.each_row() % arma::exp(0.5 * Hp.slice(n).col(t)).t();
      for (int h = 0; h <= H; ++h) out.slice((d * N + n) * (H + 1) + h) = Phi.slice(h) * P;
    }
  return out;
}

// Response of variable i to shock j at horizon h, at every date: T x N. If nv >= 0 the shock is rescaled so the
// impact response of variable nv to shock j equals `size`.
// [[Rcpp::export]]
arma::mat cpp_tvp_irf_path(const arma::cube& Bp, const arma::cube& Ap, const arma::cube& Hp, int m, int K, int p,
                           int h, int i, int j, int nv, double size) {
  const int N = Bp.n_slices, T = Bp.n_cols;
  arma::mat out(T, N);
  for (int n = 0; n < N; ++n)
    for (int t = 0; t < T; ++t) {
      arma::mat B(const_cast<double*>(Bp.slice(n).colptr(t)), m, K, false);
      arma::cube Phi = varm_ma(B, K, p, h);
      arma::mat Ainv = build_Ainv(Ap.slice(n).col(t), K);
      arma::vec shock = Ainv.col(j) * std::exp(0.5 * Hp(j, t, n));
      double sc = nv >= 0 ? size / shock(nv) : 1.0;
      out(t, n) = arma::dot(Phi.slice(h).row(i), shock) * sc;
    }
  return out;
}

// Reduced-form standard deviations sqrt(diag(Sigma_t)): K x T x N
// [[Rcpp::export]]
arma::cube cpp_tvp_sd(const arma::cube& Ap, const arma::cube& Hp, int K) {
  const int N = Ap.n_slices, T = Ap.n_cols;
  arma::cube out(K, T, N);
  for (int n = 0; n < N; ++n)
    for (int t = 0; t < T; ++t) {
      arma::mat Ainv = build_Ainv(Ap.slice(n).col(t), K);
      arma::vec lam = arma::exp(Hp.slice(n).col(t));
      for (int i = 0; i < K; ++i) out(i, t, n) = std::sqrt(arma::dot(arma::square(Ainv.row(i).t()), lam));
    }
  return out;
}

// Forecasts with simulated parameter paths. Inputs are the last-date states and hyperparameter draws of the
// stored path draws. drift = TRUE evolves beta, alpha, h forward as random walks (parameter uncertainty about the
// future); FALSE freezes them at their last values. Returns N x H x K.
// [[Rcpp::export]]
arma::cube cpp_tvp_forecast(const arma::mat& bl, const arma::mat& al, const arma::mat& hl, const arma::cube& Qp,
                            const arma::mat& Sp, const arma::mat& Wp, int m, int K, int p, const arma::mat& Ylast,
                            const arma::mat& det, int H, bool drift, bool tvb, bool tva, bool sv) {
  const int N = bl.n_cols, nb = m * K, nal = K * (K - 1) / 2, d = det.n_cols;
  arma::cube out(N, H, K);
  for (int n = 0; n < N; ++n) {
    arma::vec beta = bl.col(n), al_n = al.col(n).head(nal), h = hl.col(n);
    arma::mat Qc = chol_safe(Qp.slice(n));
    std::vector<arma::mat> Sc(K);
    int off = 0;
    for (int i = 1; i < K; ++i) {
      arma::mat Si(Sp.colptr(n) + off, i, i);
      Sc[i] = chol_safe(Si);
      off += i * i;
    }
    arma::mat Y(p + H, K);
    Y.rows(0, p - 1) = Ylast;
    arma::rowvec x(m);
    for (int s = 0; s < H; ++s) {
      if (drift) {
        if (tvb) beta += Qc * rnorm_vec(nb);
        if (tva) {
          int o = 0;
          for (int i = 1; i < K; ++i) { al_n.subvec(o, o + i - 1) += Sc[i] * rnorm_vec(i); o += i; }
        }
        if (sv) for (int i = 0; i < K; ++i) h(i) += std::sqrt(Wp(i, n)) * R::norm_rand();
      }
      arma::mat Ainv = build_Ainv(al_n, K);
      arma::vec eps(K);
      for (int i = 0; i < K; ++i) eps(i) = std::exp(0.5 * h(i)) * R::norm_rand();
      for (int l = 1; l <= p; ++l) x.subvec((l - 1) * K, l * K - 1) = Y.row(p + s - l);
      if (d > 0) x.subvec(K * p, K * p + d - 1) = det.row(s);
      arma::mat B = arma::reshape(beta, m, K);
      Y.row(p + s) = x * B + (Ainv * eps).t();
      for (int k = 0; k < K; ++k) out(n, s, k) = Y(p + s, k);
    }
  }
  return out;
}
