// [[Rcpp::depends(RcppArmadillo)]]
// Triangular ("equation-by-equation") BVAR sampler in the spirit of Carriero, Clark & Marcellino (2019):
//   A y_t = c + sum_l Pi_l y_{t-l} + Lambda_t^{1/2} eps_t,   A unit lower-triangular
//   Lambda_t = diag(exp(h_{1t}), ..., exp(h_{Kt}))       (random-walk log volatility, or constant)
// Equation i is a Gaussian regression of y_i on [x_t, -y_{1..i-1,t}] so the K equations are drawn
// independently given the volatilities: cost O(K (Kp)^3) per sweep instead of O((K^2 p)^3).
// Coefficient priors: Minnesota (fixed), SSVS (George-Sun-Ni 2008), horseshoe (Carvalho et al. 2010,
// Makalic-Schmidt 2016 representation). Volatilities: Kim-Shephard-Chib mixture (Omori et al. 2007, 10 components).
#include "varm_utils.h"
using namespace Rcpp;

namespace {

const double MIX_P[10] = {0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115};
const double MIX_M[10] = {1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000};
const double MIX_V[10] = {0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342};

inline double rig(double a, double b) { return 1.0 / R::rgamma(a, 1.0 / b); }  // inverse gamma(shape a, scale b)

// One FFBS sweep for a scalar random-walk log-volatility h (length T) given residuals e.
// Returns new h, h0 and updates sig2 (innovation variance).
void sv_update(const arma::vec& e, arma::vec& h, double& h0, double& sig2, double mu0, double V0, double a_h,
               double b_h, std::vector<int>& s) {
  const int T = e.n_elem;
  arma::vec ys(T);
  for (int t = 0; t < T; ++t) ys(t) = std::log(e(t) * e(t) + 1e-4);
  // mixture indicators
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
  // Kalman filter
  arma::vec a(T + 1), P(T + 1);
  a(0) = mu0; P(0) = V0;
  for (int t = 0; t < T; ++t) {
    double Pm = P(t) + sig2, v = MIX_V[s[t]], F = Pm + v, Kg = Pm / F;
    a(t + 1) = a(t) + Kg * (ys(t) - MIX_M[s[t]] - a(t));
    P(t + 1) = Pm * (1.0 - Kg);
  }
  // backward sampling
  h(T - 1) = a(T) + std::sqrt(P(T)) * R::norm_rand();
  for (int t = T - 2; t >= -1; --t) {
    double At = a(t + 1), Pt = P(t + 1);            // filtered moments of h_t (t=-1 -> h_0)
    double nxt = (t == T - 2) ? h(T - 1) : h(t + 1);
    double g = Pt / (Pt + sig2);
    double mean = At + g * (nxt - At), var = Pt * sig2 / (Pt + sig2);
    double draw = mean + std::sqrt(var) * R::norm_rand();
    if (t >= 0) h(t) = draw; else h0 = draw;
  }
  // innovation variance
  double ss = (h(0) - h0) * (h(0) - h0);
  for (int t = 1; t < T; ++t) ss += (h(t) - h(t - 1)) * (h(t) - h(t - 1));
  sig2 = rig(a_h + 0.5 * T, b_h + 0.5 * ss);
}

}  // namespace

// prior_type: 0 Minnesota, 1 SSVS, 2 horseshoe
// pr: [lam1, lam2, alpha, lam3 (det sd multiplier), lam5 (contemp sd), ssvs_rho, ssvs_pi]
// svp: [a_h, b_h, V0] ; homosk: [a0, b0 multiplier] prior IG(a0, b0 * s_i^2)
// [[Rcpp::export]]
List cpp_bvar_tri(const arma::mat& Y, const arma::mat& X, int p, int d, int prior_type, bool sv,
                  const arma::vec& scale, const arma::vec& delta, const arma::vec& pr, const arma::vec& svp,
                  const arma::vec& hom, int ndraw, int nburn, int thin, bool store_h) {
  const int T = Y.n_rows, K = Y.n_cols, m = X.n_cols;
  const double lam1 = pr(0), lam2 = pr(1), alpha = pr(2), lam3 = pr(3), lam5 = pr(4), rho = pr(5), pincl = pr(6);
  arma::mat Z = arma::join_rows(X, -Y);
  arma::mat G = Z.t() * Z;  // Gram of [X, -Y]: leading blocks give each equation's X'X

  // base prior variances and means (per equation, length m + i)
  std::vector<arma::vec> base(K), g0(K);
  for (int i = 0; i < K; ++i) {
    int mi = m + i;
    arma::vec bv(mi), gm(mi, arma::fill::zeros);
    for (int l = 1; l <= p; ++l)
      for (int j = 0; j < K; ++j) {
        double sc = scale(i) * scale(i) / (scale(j) * scale(j)) / std::pow((double)l, alpha);
        bv((l - 1) * K + j) = sc;
        if (l == 1 && j == i) gm((l - 1) * K + j) = delta(i);
      }
    for (int q = 0; q < d; ++q) bv(K * p + q) = std::pow(lam3 * scale(i), 2);
    for (int j = 0; j < i; ++j) bv(m + j) = std::pow(lam5, 2) * scale(i) * scale(i) / (scale(j) * scale(j));
    base[i] = bv; g0[i] = gm;
  }
  const int nlag = K * p;  // shrinkable block is indices 0..nlag-1
  // state
  arma::mat Pi(m, K), A(K, K, arma::fill::eye);
  arma::mat H(T, K);
  for (int i = 0; i < K; ++i) H.col(i).fill(std::log(scale(i) * scale(i)));
  arma::vec h0(K), sig2h(K), sig2e(K);
  for (int i = 0; i < K; ++i) { h0(i) = std::log(scale(i) * scale(i)); sig2h(i) = svp(1) / (svp(0) - 1); sig2e(i) = scale(i) * scale(i); }
  // shrinkage state
  std::vector<arma::vec> lam2v(K), nuv(K), delt(K);
  arma::vec tau2(K, arma::fill::ones), xi(K, arma::fill::ones);
  for (int i = 0; i < K; ++i) { lam2v[i].ones(nlag); nuv[i].ones(nlag); delt[i].ones(nlag); tau2(i) = lam1 * lam1; }
  std::vector<int> s(T);

  const int nsave = ndraw;
  arma::cube Bs(m, K, nsave), As(K, K, nsave), Sl(K, K, nsave), Sm(K, K, nsave);
  arma::mat hlast(nsave, K), sh(nsave, K), hmean(T, K, arma::fill::zeros), incl(nlag, K, arma::fill::zeros);
  arma::mat tauS(nsave, K);
  arma::cube hall;
  if (store_h) hall.set_size(nsave, T, K);
  int saved = 0;
  const int total = nburn + ndraw * thin;
  for (int it = 0; it < total; ++it) {
    for (int i = 0; i < K; ++i) {
      const int mi = m + i;
      // prior variance vector
      arma::vec v = base[i];
      switch (prior_type) {
        case 0:
          for (int j = 0; j < nlag; ++j) v(j) = base[i](j) * lam1 * lam1 * (j % K == i ? 1.0 : lam2 * lam2);
          break;
        case 1:
          for (int j = 0; j < nlag; ++j) {
            double v1 = base[i](j) * lam1 * lam1;
            v(j) = delt[i](j) > 0.5 ? v1 : rho * rho * v1;
          }
          break;
        default:
          for (int j = 0; j < nlag; ++j) v(j) = base[i](j) * tau2(i) * lam2v[i](j);
      }
      // posterior of gamma_i = (pi_i, a_i,1..i-1)
      arma::mat XtX; arma::vec Xty;
      if (sv) {
        arma::vec w = arma::exp(-0.5 * H.col(i));
        arma::mat Xi = Z.cols(0, mi - 1);
        arma::mat Xw = Xi.each_col() % w;
        XtX = Xw.t() * Xw;
        Xty = Xw.t() * (Y.col(i) % w);
      } else {
        XtX = G.submat(0, 0, mi - 1, mi - 1) / sig2e(i);
        Xty = -G.submat(0, m + i, mi - 1, m + i) / sig2e(i);
      }
      arma::vec prec = 1.0 / v;
      arma::mat Pm = XtX; Pm.diag() += prec;
      arma::mat Lc = arma::chol(Pm, "lower");
      arma::vec rhs = Xty + prec % g0[i];
      arma::vec mu = arma::solve(arma::trimatu(Lc.t()), arma::solve(arma::trimatl(Lc), rhs));
      arma::vec zz(mi);
      zz.imbue([]() { return R::norm_rand(); });
      arma::vec gam = mu + arma::solve(arma::trimatu(Lc.t()), zz);
      Pi.col(i) = gam.head(m);
      for (int j = 0; j < i; ++j) A(i, j) = gam(m + j);
      // residuals
      arma::vec e = Y.col(i) - Z.cols(0, mi - 1) * gam;
      // shrinkage updates (lag block)
      if (prior_type == 1) {
        for (int j = 0; j < nlag; ++j) {
          double v1 = base[i](j) * lam1 * lam1, v0 = rho * rho * v1, x = gam(j) - g0[i](j);
          double l1 = std::log(pincl) - 0.5 * std::log(v1) - 0.5 * x * x / v1;
          double l0 = std::log(1 - pincl) - 0.5 * std::log(v0) - 0.5 * x * x / v0;
          double pr1 = 1.0 / (1.0 + std::exp(l0 - l1));
          delt[i](j) = R::unif_rand() < pr1 ? 1.0 : 0.0;
        }
      } else if (prior_type == 2) {
        double sumt = 0;
        for (int j = 0; j < nlag; ++j) {
          double x = gam(j) - g0[i](j), th = x * x / base[i](j);
          lam2v[i](j) = rig(1.0, 1.0 / nuv[i](j) + th / (2 * tau2(i)));
          nuv[i](j) = rig(1.0, 1.0 + 1.0 / lam2v[i](j));
          sumt += th / lam2v[i](j);
        }
        tau2(i) = rig(0.5 * (nlag + 1), 1.0 / xi(i) + 0.5 * sumt);
        xi(i) = rig(1.0, 1.0 + 1.0 / tau2(i));
      }
      // variance / volatility update
      if (sv) {
        arma::vec h = H.col(i);
        sv_update(e, h, h0(i), sig2h(i), std::log(scale(i) * scale(i)), svp(2), svp(0), svp(1), s);
        H.col(i) = h;
      } else {
        sig2e(i) = rig(hom(0) + 0.5 * T, hom(1) * scale(i) * scale(i) + 0.5 * arma::dot(e, e));
        H.col(i).fill(std::log(sig2e(i)));
      }
    }
    if (it >= nburn && ((it - nburn) % thin == thin - 1)) {
      arma::mat Ainv = arma::inv(arma::trimatl(A));
      Bs.slice(saved) = Pi * Ainv.t();
      As.slice(saved) = A;
      arma::vec hl = H.row(T - 1).t();
      Sl.slice(saved) = Ainv * arma::diagmat(arma::exp(hl)) * Ainv.t();
      arma::vec ev = arma::mean(arma::exp(H), 0).t();
      Sm.slice(saved) = Ainv * arma::diagmat(ev) * Ainv.t();
      hlast.row(saved) = hl.t();
      sh.row(saved) = arma::sqrt(sig2h).t();
      tauS.row(saved) = arma::sqrt(tau2).t();
      hmean += H;
      if (store_h) for (int k = 0; k < K; ++k) hall.slice(k).row(saved) = H.col(k).t();
      if (prior_type == 1) for (int i = 0; i < K; ++i) incl.col(i) += delt[i];
      ++saved;
    }
    if ((it & 63) == 0) Rcpp::checkUserInterrupt();
  }
  hmean /= nsave; incl /= nsave;
  List out = List::create(_["B"] = Bs, _["A"] = As, _["Sigma"] = Sl, _["Sigma_mean"] = Sm,
                          _["h_last"] = hlast, _["sigma_h"] = sh, _["h_mean"] = hmean,
                          _["inclusion"] = incl, _["tau"] = tauS);
  if (store_h) out["h_all"] = hall;
  return out;
}

// Forecast with stochastic-volatility simulation. For each posterior draw, simulate h forward
// as a random walk with innovation sd sigma_h, reduced-form shocks u = A^{-1} Lambda^{1/2} eps.
// Returns cube [draw, h, var].
// [[Rcpp::export]]
arma::cube cpp_forecast_sv(const arma::cube& Bs, const arma::cube& As, const arma::mat& hlast,
                           const arma::mat& sigh, int p, const arma::mat& Ylast, const arma::mat& det, int H) {
  const int N = Bs.n_slices, K = Ylast.n_cols;
  arma::cube out(N, H, K);
  for (int n = 0; n < N; ++n) {
    arma::mat Ainv = arma::inv(arma::trimatl(As.slice(n)));
    arma::vec h = hlast.row(n).t();
    arma::mat U(H, K);
    for (int s = 0; s < H; ++s) {
      arma::vec eps(K);
      for (int k = 0; k < K; ++k) {
        h(k) += sigh(n, k) * R::norm_rand();
        eps(k) = std::exp(0.5 * h(k)) * R::norm_rand();
      }
      U.row(s) = (Ainv * eps).t();
    }
    arma::mat path = varm_path(Bs.slice(n), K, p, Ylast, det, U);
    for (int k = 0; k < K; ++k) out.slice(k).row(n) = path.col(k).t();
  }
  return out;
}
