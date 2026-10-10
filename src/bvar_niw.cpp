// [[Rcpp::depends(RcppArmadillo)]]
// Conjugate Normal-inverse-Wishart BVAR with Minnesota prior, sum-of-coefficients and
// single-unit-root (dummy initial observation) priors, closed-form marginal likelihood and
// hierarchical hyperparameter sampling (Giannone, Lenza & Primiceri 2015).
#include "varm_utils.h"
using namespace Rcpp;

namespace {

struct NIWSpec {
  int M, p, d;        // variables, lags, deterministic columns
  bool has_const;
  arma::mat Y, X;     // regression sample (T x M, T x K)
  arma::vec ybar;     // mean of the first p observations (for dummy priors)
  arma::mat b;        // prior mean (K x M)
  double vard;        // prior variance for deterministic coefficients
};

struct NIWPost {
  double logml;
  arma::mat bhat, S, Ginv;  // posterior mean, IW scale (psi+sse+penalty), (X'X+Om^-1)^-1
  double nu;                // IW degrees of freedom (N + M + 2)
};

// Evaluate marginal likelihood (and posterior) for hyperparameters
// lambda, alpha, psi (M), mu (<=0: off), th (<=0: off).
NIWPost niw_eval(const NIWSpec& s, double lambda, double alpha, const arma::vec& psi, double mu, double th,
                 bool full = true) {
  const int M = s.M, p = s.p, K = M * p + s.d;
  arma::vec om(K);
  for (int l = 1; l <= p; ++l)
    for (int i = 0; i < M; ++i) om((l - 1) * M + i) = lambda * lambda / std::pow((double)l, alpha) / psi(i);
  for (int j = 0; j < s.d; ++j) om(M * p + j) = s.vard;

  // dummy observations
  arma::mat Yd(0, M), Xd(0, K);
  if (mu > 0) {
    arma::mat soc = arma::diagmat(s.ybar) / mu;
    arma::mat Xs(M, K, arma::fill::zeros);
    for (int l = 0; l < p; ++l) Xs.cols(l * M, (l + 1) * M - 1) = soc;
    Yd = arma::join_cols(Yd, soc);
    Xd = arma::join_cols(Xd, Xs);
  }
  if (th > 0 && s.has_const) {
    arma::rowvec ys = s.ybar.t() / th;
    arma::rowvec xs(K, arma::fill::zeros);
    for (int l = 0; l < p; ++l) xs.subvec(l * M, (l + 1) * M - 1) = ys;
    xs(M * p) = 1.0 / th;   // constant is the first deterministic column
    Yd = arma::join_cols(Yd, ys);
    Xd = arma::join_cols(Xd, xs);
  }
  const int Nd = Yd.n_rows;
  arma::mat Ya = arma::join_cols(Yd, s.Y), Xa = arma::join_cols(Xd, s.X);
  const int N = Ya.n_rows;
  arma::mat Oinv = arma::diagmat(1.0 / om), Osq = arma::diagmat(arma::sqrt(om));
  arma::vec spsi = arma::sqrt(psi);
  arma::mat psi_inv_sq = arma::diagmat(1.0 / spsi);

  auto ev = [&](const arma::mat& Yv, const arma::mat& Xv, bool est, arma::mat& bh, arma::mat& sse,
                arma::mat& pen, arma::vec& oev, arma::vec& pev) {
    arma::mat XX = Xv.t() * Xv;
    if (est) bh = arma::solve(XX + Oinv, Xv.t() * Yv + Oinv * s.b, arma::solve_opts::likely_sympd);
    else bh = s.b;
    arma::mat R = Yv - Xv * bh;
    sse = R.t() * R;
    pen = est ? arma::mat((bh - s.b).t() * Oinv * (bh - s.b)) : arma::mat(M, M, arma::fill::zeros);
    arma::mat oml = Osq * XX * Osq;
    arma::mat pml = psi_inv_sq * (sse + pen) * psi_inv_sq;
    oev = arma::eig_sym(oml);
    pev = arma::eig_sym(pml);
    oev.transform([](double x) { return (x < 1e-12 ? 0.0 : x) + 1.0; });
    pev.transform([](double x) { return (x < 1e-12 ? 0.0 : x) + 1.0; });
  };
  auto lml = [&](int Nn, const arma::vec& oev, const arma::vec& pev) {
    double r = -0.5 * M * Nn * std::log(M_PI);
    for (int j = 0; j < M; ++j) r += R::lgammafn((Nn + M + 2 - j) / 2.0) - R::lgammafn((M + 2 - j) / 2.0);
    r -= 0.5 * Nn * arma::accu(arma::log(psi));
    r -= 0.5 * M * arma::accu(arma::log(oev));
    r -= 0.5 * (Nn + M + 2) * arma::accu(arma::log(pev));
    return r;
  };

  NIWPost out;
  arma::mat bh, sse, pen; arma::vec oev, pev;
  ev(Ya, Xa, true, bh, sse, pen, oev, pev);
  out.logml = lml(N, oev, pev);
  if (Nd > 0) {
    arma::mat b2, sse2, pen2; arma::vec oev2, pev2;
    ev(Yd, Xd, false, b2, sse2, pen2, oev2, pev2);
    out.logml -= lml(Nd, oev2, pev2);
  }
  out.nu = N + M + 2;
  if (!full) return out;   // posterior moments are only needed for accepted / retained points
  out.bhat = bh;
  out.S = arma::diagmat(psi) + sse + pen;
  out.Ginv = arma::inv_sympd(Xa.t() * Xa + Oinv);
  return out;
}

// Draw (B, Sigma) ~ NIW posterior.
void niw_draw(const NIWPost& P, arma::mat& B, arma::mat& Sig) {
  const int M = P.S.n_rows, K = P.bhat.n_rows;
  // Sigma ~ IW(S, nu): Sigma = (L_Sinv A A' L_Sinv')^-1 with Bartlett A
  arma::mat LSi = arma::chol(arma::inv_sympd(P.S), "lower");
  arma::mat A(M, M, arma::fill::zeros);
  for (int i = 0; i < M; ++i) {
    A(i, i) = std::sqrt(R::rchisq(P.nu - i));
    for (int j = 0; j < i; ++j) A(i, j) = R::norm_rand();
  }
  arma::mat LW = LSi * A;
  Sig = arma::inv_sympd(LW * LW.t());
  arma::mat LS = arma::chol(Sig, "lower");
  arma::mat LG = arma::chol(P.Ginv, "lower");
  arma::mat Z(K, M);
  Z.imbue([]() { return R::norm_rand(); });
  B = P.bhat + LG * Z * LS.t();
}

NIWSpec make_spec(const arma::mat& Y, const arma::mat& X, int p, int d, bool has_const, const arma::mat& b,
                  double vard) {
  NIWSpec s;
  s.M = Y.n_cols; s.p = p; s.d = d; s.has_const = has_const;
  s.Y = Y; s.X = X; s.b = b; s.vard = vard;
  s.ybar = arma::mean(Y.rows(0, p - 1), 0).t();
  return s;
}

double gamma_logpdf(double x, double shape, double scale) { return R::dgamma(x, shape, scale, 1); }
double ig_logpdf(double x, double shape, double scale) {
  return shape * std::log(scale) - (shape + 1) * std::log(x) - scale / x - R::lgammafn(shape);
}

// hyper layout: [lambda, alpha, mu, theta, psi_1..psi_M]
double log_post(const NIWSpec& s, const arma::vec& h, const arma::uvec& isfree, const arma::mat& gp,
                const arma::vec& lo, const arma::vec& hi, const arma::vec& ig, bool eval_only_ml, NIWPost* post) {
  for (arma::uword i = 0; i < h.n_elem; ++i)
    if (h(i) < lo(i) || h(i) > hi(i)) return -1e18;
  arma::vec psi = h.subvec(4, h.n_elem - 1);
  NIWPost P = niw_eval(s, h(0), h(1), psi, h(2), h(3), post != nullptr);
  if (post) *post = P;
  double lp = P.logml;
  if (eval_only_ml) return lp;
  for (int i = 0; i < 4; ++i)
    if (isfree(i)) lp += gamma_logpdf(h(i), gp(i, 0), gp(i, 1));
  for (arma::uword i = 4; i < h.n_elem; ++i)
    if (isfree(i)) lp += ig_logpdf(h(i), ig(0), ig(1));
  return lp;
}

}  // namespace

// Log marginal likelihood (+ hyperprior if eval_only_ml = false) at a given hyperparameter vector
// [[Rcpp::export]]
double cpp_niw_logpost(const arma::vec& h, const arma::mat& Y, const arma::mat& X, int p, int d,
                       bool has_const, const arma::mat& b, double vard, const arma::uvec& isfree,
                       const arma::mat& gp, const arma::vec& lo, const arma::vec& hi, const arma::vec& ig,
                       bool eval_only_ml) {
  NIWSpec s = make_spec(Y, X, p, d, has_const, b, vard);
  return log_post(s, h, isfree, gp, lo, hi, ig, eval_only_ml, nullptr);
}

// Posterior draws at fixed hyperparameters (direct sampling: no MCMC needed)
// [[Rcpp::export]]
List cpp_niw_fixed(const arma::vec& h, const arma::mat& Y, const arma::mat& X, int p, int d, bool has_const,
                   const arma::mat& b, double vard, int ndraw) {
  NIWSpec s = make_spec(Y, X, p, d, has_const, b, vard);
  arma::vec psi = h.subvec(4, h.n_elem - 1);
  NIWPost P = niw_eval(s, h(0), h(1), psi, h(2), h(3));
  const int K = X.n_cols, M = Y.n_cols;
  arma::cube Bs(K, M, ndraw), Ss(M, M, ndraw);
  arma::mat B, S;
  for (int n = 0; n < ndraw; ++n) {
    niw_draw(P, B, S);
    Bs.slice(n) = B; Ss.slice(n) = S;
  }
  return List::create(_["B"] = Bs, _["Sigma"] = Ss, _["logml"] = P.logml, _["bhat"] = P.bhat, _["S"] = P.S,
                      _["nu"] = P.nu);
}

// Hierarchical BVAR: random-walk Metropolis on the free hyperparameters, one (B, Sigma) draw per
// retained iteration. Proposal covariance `prop` (nfree x nfree) is scaled adaptively during burn-in.
// [[Rcpp::export]]
List cpp_niw_hier(const arma::vec& h0, const arma::mat& Y, const arma::mat& X, int p, int d, bool has_const,
                  const arma::mat& b, double vard, const arma::uvec& isfree, const arma::mat& gp,
                  const arma::vec& lo, const arma::vec& hi, const arma::vec& ig, const arma::mat& prop,
                  int niter, int nburn, double scale, bool adapt, double acc_lo, double acc_hi) {
  NIWSpec s = make_spec(Y, X, p, d, has_const, b, vard);
  const arma::uvec idx = arma::find(isfree);
  const int nf = idx.n_elem, K = X.n_cols, M = Y.n_cols;
  arma::mat Lp = arma::chol(prop, "lower");
  arma::vec h = h0;
  NIWPost P;
  double lp = log_post(s, h, isfree, gp, lo, hi, ig, false, &P);
  arma::cube Bs(K, M, niter), Ss(M, M, niter);
  arma::mat Hs(niter, h.n_elem), B, S;
  arma::vec lps(niter);
  long acc = 0;
  int blk = 0; long blk_acc = 0;
  for (int it = 0; it < niter + nburn; ++it) {
    arma::vec z(nf);
    z.imbue([]() { return R::norm_rand(); });
    arma::vec hn = h;
    hn.elem(idx) += scale * (Lp * z);
    double lpn = log_post(s, hn, isfree, gp, lo, hi, ig, false, nullptr);
    bool ok = std::log(R::unif_rand()) < lpn - lp;
    if (ok) { h = hn; lp = lpn; log_post(s, hn, isfree, gp, lo, hi, ig, false, &P); ++blk_acc; if (it >= nburn) ++acc; }
    ++blk;
    if (adapt && it < nburn && blk == 50) {
      double r = (double)blk_acc / blk;
      if (r < acc_lo) scale *= 0.9; else if (r > acc_hi) scale *= 1.1;
      blk = 0; blk_acc = 0;
    }
    if (it >= nburn) {
      int k = it - nburn;
      niw_draw(P, B, S);
      Bs.slice(k) = B; Ss.slice(k) = S;
      Hs.row(k) = h.t();
      lps(k) = lp;
    }
    if ((it & 255) == 0) Rcpp::checkUserInterrupt();
  }
  return List::create(_["B"] = Bs, _["Sigma"] = Ss, _["hyper"] = Hs, _["logpost"] = lps,
                      _["accept"] = (double)acc / niter, _["scale"] = scale);
}
