// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <R_ext/Utils.h>
#include <R_ext/Print.h>

using namespace Rcpp;
using namespace arma;

inline bool safe_chol(mat& L, const mat& A, double& ridge, const bool lower=true) {
  const double diag_mean = mean(A.diag());
  ridge = std::max(ridge, 1e-12 * std::max(1.0, std::abs(diag_mean)));
  mat M;
  for (int i = 0; i < 8; ++i) {
    M = A + ridge * eye<mat>(A.n_rows, A.n_cols);
    if (chol(L, M, lower ? "lower" : "upper")) return true;
    ridge *= 10.0;
  }
  return false;
}

inline vec spd_solve_from_chol(const mat& L, const vec& s) {
  vec v = solve(trimatl(L),     s, solve_opts::fast);
  vec z = solve(trimatu(L.t()), v, solve_opts::fast);
  return z;
}

static double log_post_SL(const vec& beta,
                          const vec& b0,
                          const mat& B_inv,
                          const vec& y,
                          const mat& X,
                          const vec& w,
                          double tau,
                          const mat& L_XtW2X) {
  vec diff = beta - b0;
  double lp = -0.5 * dot(diff, B_inv * diff);

  vec res = y - X * beta;
  vec ind = tau - conv_to<vec>::from(res < 0);
  vec s   = X.t() * (w % ind);

  vec z = spd_solve_from_chol(L_XtW2X, s);
  double quad = 0.5 / (tau * (1.0 - tau)) * dot(s, z);

  return lp - quad;
}

Rcpp::List _mcmc_bwqr_sl_cpp(const arma::vec& y,
                             const arma::mat& X,
                             const arma::vec& w,
                             double tau = 0.5,
                             int n_mcmc  = 10000,
                             int burnin  = 2000,
                             int thin    = 10,
                             Rcpp::Nullable<Rcpp::NumericVector> b_prior_mean = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericMatrix> B_prior_prec = R_NilValue,
                             int print_progress = 0) {
  if (y.n_elem != X.n_rows || w.n_elem != y.n_elem)
    stop("Incompatible dimensions between y, X, and w.");
  if (!(tau > 0.0 && tau < 1.0))
    stop("tau must be in (0,1).");
  if (burnin >= n_mcmc) stop("burnin must be < n_mcmc.");
  if (thin <= 0)        stop("thin must be positive.");

  const int n = y.n_elem;
  const int p = X.n_cols;

  vec b0 = arma::zeros<vec>(p);
  if (b_prior_mean.isNotNull()) {
    b0 = as<vec>(b_prior_mean);
    if ((int)b0.n_elem != p) stop("b_prior_mean must have length equal to ncol(X).");
  }

  mat B_inv = eye<mat>(p, p) / 1000.0;
  if (B_prior_prec.isNotNull()) {
    B_inv = as<mat>(B_prior_prec);
    if (B_inv.n_rows != (uword)p || B_inv.n_cols != (uword)p)
      stop("B_prior_prec must be a p x p matrix.");
  }

  if (w.min() <= 0.0) stop("All weights must be > 0.");

  vec w2 = w % w;
  mat XtW2X = X.t() * (X.each_col() % w2);
  mat L_XtW2X;
  double ridge_X = 0.0;
  if (!safe_chol(L_XtW2X, XtW2X, ridge_X, true))
    stop("Cholesky decomposition of X' W^2 X failed even with ridge.");

  mat I_p = eye<mat>(p, p);
  mat L_prop = solve(trimatu(L_XtW2X.t()), I_p, solve_opts::fast);
  const double prop_scalar = std::sqrt(tau * (1.0 - tau) * (double)n);

  const int n_keep = (n_mcmc - burnin) / thin;
  mat beta_out(n_keep, p, fill::none);
  int accept = 0, k_out = 0;

  vec beta = solve(X, y);
  double ct = 1.0;

  double logp_curr = log_post_SL(beta, b0, B_inv, y, X, w, tau, L_XtW2X);

  const int bar_width = 40;
  int last_decile = -1;

  for (int k = 0; k < n_mcmc; ++k) {
    if (print_progress > 0) {
      double perc = 100.0 * (k + 1.0) / n_mcmc;
      int decile = (static_cast<int>(std::floor(perc)) / 10) * 10;
      if (decile >= 10 && decile <= 100 && decile > last_decile) {
        int filled = static_cast<int>(std::round(bar_width * perc / 100.0));
        Rprintf("\r[");
        for (int j = 0; j < bar_width; ++j) Rprintf(j < filled ? "=" : " ");
        Rprintf("] %5.1f%%", perc);
        R_FlushConsole();
        R_CheckUserInterrupt();
        last_decile = decile;
      }
    }

    vec z = randn<vec>(p);
    vec beta_prop = beta + std::sqrt(ct) * prop_scalar * (L_prop * z);

    double logp_prop = log_post_SL(beta_prop, b0, B_inv, y, X, w, tau, L_XtW2X);

    double log_alpha = std::min(0.0, logp_prop - logp_curr);
    double u = R::runif(0.0, 1.0);
    const double acc_prob = std::exp(log_alpha);

    if (u < acc_prob) {
      beta = beta_prop;
      logp_curr = logp_prop;
      ++accept;
    }

    double step = std::pow(k + 1.0, -0.8);
    ct = std::exp(std::log(ct) + step * (acc_prob - 0.234));
    ct = std::max(1e-4, std::min(ct, 1e4));

    if (k >= burnin && ((k - burnin) % thin == 0)) {
      beta_out.row(k_out++) = beta.t();
    }
  }

  if (print_progress > 0) {
    Rprintf("\r[");
    for (int j = 0; j < bar_width; ++j) Rprintf("=");
    Rprintf("] %5.1f%%\n", 100.0);
    R_FlushConsole();
  }

  return List::create(
    _["beta"]        = beta_out,
    _["accept_rate"] = double(accept) / n_mcmc,
    _["call"]        = "MCMC_BWQR_SL"
  );
}
