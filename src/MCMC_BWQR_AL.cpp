// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <R_ext/Utils.h>   // R_CheckUserInterrupt
#include <R_ext/Print.h>   // Rprintf, R_FlushConsole

using namespace Rcpp;
using namespace arma;

inline double rinvgauss(double mu, double lambda) {
  const double z = R::rnorm(0.0, 1.0);
  const double y = z * z;
  const double x1 = mu + (mu * mu * y) / (2.0 * lambda) -
    (mu / (2.0 * lambda)) * std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
  const double u = R::runif(0.0, 1.0);
  return (u <= mu / (mu + x1)) ? x1 : (mu * mu / x1);
}


static arma::vec draw_beta(
    const arma::mat& X,
    const arma::vec& w,
    const arma::vec& v,
    const arma::vec& y,
    double sigma,
    double delta2,
    double theta,
    const arma::mat& B_inv,
    const arma::vec& b_mean) {

  const double sigma_eff = std::max(sigma, 1e-12);
  const arma::vec wv = w / v;
  arma::mat XtWvX = X.t() * (X.each_col() % wv);
  arma::mat A = symmatu(B_inv + XtWvX / (delta2 * sigma_eff));

  arma::mat Sigma;
  double ridge = 1e-8;
  bool ok = false;
  for (int tries = 0; tries < 6 && !ok; ++tries) {
    ok = arma::inv_sympd(Sigma, A + ridge * arma::eye<arma::mat>(A.n_rows, A.n_cols));
    if (!ok) ridge *= 10.0;
  }
  if (!ok) {
    Sigma = arma::pinv(A, 1e-12);
    ridge = 1e-6;
  }

  // prior term
  arma::vec rhs = B_inv * b_mean + X.t() * (wv % (y - theta * v)) / (delta2 * sigma_eff);
  arma::vec mu = Sigma * rhs;

  arma::mat L;
  ok = arma::chol(L, Sigma + ridge * arma::eye<arma::mat>(Sigma.n_rows, Sigma.n_cols), "lower");
  if (!ok) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Sigma);
    arma::mat Sigma_pd = eigvec * arma::diagmat(arma::clamp(eigval, 1e-12, arma::datum::inf)) * eigvec.t();
    arma::chol(L, Sigma_pd, "lower");
  }

  return mu + L * arma::randn<arma::vec>(X.n_cols);
}


inline double draw_sigma(
    const arma::mat& X,
    const arma::vec& w,
    const arma::vec& beta,
    const arma::vec& v,
    const arma::vec& y,
    double tau2,
    double theta,
    double c0,
    double C0) {

  const int n = y.n_elem;
  const double alpha1 = c0 + 1.5 * n;
  arma::vec resid = y - X * beta - theta * v;
  const double quad = arma::dot((w % resid) / v, resid) / (2.0 * tau2);
  const double beta1 = C0 + arma::sum(w % v) + quad;
  // R::rgamma(shape, scale); scale = 1/rate
  return 1.0 / R::rgamma(alpha1, 1.0 / beta1);
}


static void update_v(
    arma::vec& v,
    const arma::mat& X,
    const arma::vec& w,
    const arma::vec& beta,
    const arma::vec& y,
    double delta2,
    double theta,
    double sigma) {

  const int N = v.n_elem;
  const double sigma_eff = std::max(sigma, 1e-12);

  for (int i = 0; i < N; ++i) {
    const double resid = y[i] - arma::dot(X.row(i), beta);
    const double chi = w[i] * resid * resid / (delta2 * sigma_eff);
    const double psi = w[i] * (2.0 / sigma_eff + theta * theta / (delta2 * sigma_eff));
    const double mu = std::sqrt(chi / psi);
    const double lambda = std::sqrt(chi * psi);
    v[i] = std::max(rinvgauss(mu, lambda), 1e-8);
  }
}

Rcpp::List _mcmc_bwqr_al_cpp(
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& w,
    double tau = 0.5,
    int n_mcmc = 50000,
    int burnin = 10000,
    int thin = 10,
    Rcpp::Nullable<Rcpp::NumericVector> b_prior_mean = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> B_prior_prec = R_NilValue,
    double c0 = 0.001,
    double C0 = 0.001,
    int print_progress = 1000,
    Rcpp::Nullable<Rcpp::NumericVector> fix_sigma = R_NilValue
) {

  if (y.n_elem != X.n_rows || y.n_elem != w.n_elem)
    stop("Dimensions of y, X, and w must match.");
  if (burnin >= n_mcmc)
    stop("burnin must be < n_mcmc.");
  if (thin <= 0)
    stop("thin must be positive.");

  const int n = y.n_elem, p = X.n_cols;
  arma::vec b_mean = b_prior_mean.isNotNull() ? Rcpp::as<arma::vec>(b_prior_mean) : arma::zeros<arma::vec>(p);
  if ((int)b_mean.n_elem != p)
    stop("b_prior_mean must have length equal to ncol(X).");
  arma::mat B_prec = B_prior_prec.isNotNull() ?
    Rcpp::as<arma::mat>(B_prior_prec) : (arma::eye<arma::mat>(p, p) / 1000.0);
  if ((int)B_prec.n_rows != p || (int)B_prec.n_cols != p)
    stop("B_prior_prec must be a p x p matrix.");

  bool use_fixed_sigma = false;
  double sigma_fixed_val = 1.0;
  if (fix_sigma.isNotNull()) {
    NumericVector fs(fix_sigma);
    if (fs.size() != 1 || !R_finite(fs[0]) || fs[0] <= 0.0)
      stop("'fix_sigma' must be a positive scalar if provided.");
    use_fixed_sigma = true;
    sigma_fixed_val = fs[0];
  }

  arma::mat beta_chain(n_mcmc, p, arma::fill::zeros);
  arma::vec sigma_chain(n_mcmc, arma::fill::zeros);

  beta_chain.row(0) = arma::solve(X, y).t();
  sigma_chain[0] = use_fixed_sigma ? sigma_fixed_val : 1.0;
  arma::vec v = arma::randg<arma::vec>(n, arma::distr_param(2.0, 1.0));

  const double delta2 = 2.0 / (tau * (1.0 - tau));
  const double theta  = (1.0 - 2.0 * tau) / (tau * (1.0 - tau));
  const double tau2   = delta2;

  const int bar_width = 40;
  int last_decile = -1;

  for (int k = 1; k < n_mcmc; ++k) {
    if (print_progress > 0) {
      double perc = 100.0 * (k + 0.0) / (n_mcmc - 1.0);
      int decile = ((int)std::floor(perc)) / 10 * 10;
      if (decile >= 10 && decile <= 100 && decile > last_decile) {
        int filled = (int)std::round(bar_width * perc / 100.0);
        Rprintf("\r[");
        for (int j = 0; j < bar_width; ++j) Rprintf(j < filled ? "=" : " ");
        Rprintf("] %5.1f%%", perc);
        R_FlushConsole();
        R_CheckUserInterrupt();
        last_decile = decile;
      }
    }

    arma::vec beta_k = draw_beta(X, w, v, y, sigma_chain[k - 1], delta2, theta, B_prec, b_mean);
    beta_chain.row(k) = beta_k.t();

    update_v(v, X, w, beta_k, y, delta2, theta, sigma_chain[k - 1]);

    if (use_fixed_sigma) {
      sigma_chain[k] = sigma_fixed_val;
    } else {
      sigma_chain[k] = draw_sigma(X, w, beta_k, v, y, tau2, theta, c0, C0);
    }
  }

  if (print_progress > 0) {
    Rprintf("\r[");
    for (int j = 0; j < bar_width; ++j) Rprintf("=");
    Rprintf("] %5.1f%%\n", 100.0);
    R_FlushConsole();
  }

  arma::uvec keep = arma::regspace<arma::uvec>(burnin, thin, n_mcmc - 1);
  const int M = keep.n_elem;

  NumericMatrix beta_out(M, p);
  NumericVector sigma_out(M);

  for (int i = 0; i < M; ++i) {
    const int ki = keep[i];
    for (int j = 0; j < p; ++j) beta_out(i, j) = beta_chain(ki, j);
    sigma_out[i] = sigma_chain[ki];
  }

  return List::create(
    _["beta"]  = beta_out,
    _["sigma"] = sigma_out
  );
}
