// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <R_ext/Utils.h>   // R_CheckUserInterrupt
#include <R_ext/Print.h>   // Rprintf, R_FlushConsole

using namespace Rcpp;
using namespace arma;

constexpr double PI = 3.14159265358979323846;

// --- Safe Cholesky with growing ridge ---
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
  vec v = solve(trimatl(L), s, solve_opts::fast);     // L v = s
  vec z = solve(trimatl(L.t()), v, solve_opts::fast); // L' z = v
  return z;
}

static double log_post_SL_stable(const vec& beta,
                                 const vec& b0,
                                 const mat& B_inv,
                                 const vec& y,
                                 const mat& X,
                                 const vec& w,        // already normalized in R
                                 double tau,
                                 const mat& L_XtWX,   // chol (lower) of XtWX (with normalized w)
                                 double lambda) {
  // Prior
  vec diff = beta - b0;
  double lp = -0.5 * dot(diff, B_inv * diff);

  // Score s_tau
  vec res  = y - X * beta;
  vec ind  = tau - conv_to<vec>::from(res < 0);
  vec s    = X.t() * (w % ind);

  vec z = spd_solve_from_chol(L_XtWX, s);
  double quad = 0.5 * (lambda / (tau * (1.0 - tau))) * dot(s, z);

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
  // --- Basic checks ---
  if (y.n_elem != X.n_rows || w.n_elem != y.n_elem)
    stop("Incompatible dimensions between y, X, and w.");
  if (!(tau > 0.0 && tau < 1.0))
    stop("tau must be in (0,1).");
  if (burnin >= n_mcmc) stop("burnin must be < n_mcmc.");
  if (thin <= 0)        stop("thin must be positive.");

  const int p = X.n_cols;

  // --- Prior ---
  vec b0 = arma::zeros<vec>(p);
  if (b_prior_mean.isNotNull()) {
    b0 = as<vec>(b_prior_mean);
    if ((int)b0.n_elem != p) stop("b_prior_mean must have length equal to ncol(X).");
  }

  mat B_inv = eye<mat>(p, p) / 1e6; // very diffuse prior by default
  if (B_prior_prec.isNotNull()) {
    B_inv = as<mat>(B_prior_prec);
    if (B_inv.n_rows != (uword)p || B_inv.n_cols != (uword)p)
      stop("B_prior_prec must be a p x p matrix.");
  }

  // --- Weights: use as they come (already normalized in R) ---
  if (w.min() <= 0.0) stop("All weights must be > 0.");
  const double mean_w = arma::mean(w);
  if (std::abs(mean_w - 1.0) > 1e-8) {
    Rcpp::warning("Expected weights normalized to mean 1; mean(w)=%.6f", mean_w);
  }
  const double lambda = sum(w); // ≈ n if mean(w)=1

  // --- XtWX and its safe Cholesky ---
  mat XtWX = X.t() * (X.each_col() % w);
  mat L_XtWX;
  double ridge_X = 0.0;
  if (!safe_chol(L_XtWX, XtWX, ridge_X, true))
    stop("Cholesky decomposition of XtWX failed even with ridge.");

  // --- Proposal: preconditioned and stable ---
  const double kappa = 1.0;
  mat Prec = B_inv + kappa * XtWX;
  mat L_Prec;
  double ridge_P = 0.0;
  if (!safe_chol(L_Prec, Prec, ridge_P, true))
    stop("Cholesky decomposition of (B_inv + kappa*XtWX) failed.");

  mat I_p = eye<mat>(p, p);
  mat L_prop_base = solve(trimatl(L_Prec.t()), I_p, solve_opts::fast);
  double scale0 = (p > 1) ? (2.38 * 2.38 / p) : 1.0;

  // --- Output ---
  const int n_keep = (n_mcmc - burnin) / thin;
  mat beta_out(n_keep, p, fill::none);
  int accept = 0, k_out = 0;

  // --- Init ---
  vec beta = solve(X, y); // LS
  double ct = 1.0;

  // Precompute current log-posterior
  double logp_curr = log_post_SL_stable(beta, b0, B_inv, y, X, w, tau, L_XtWX, lambda);

  const int bar_width = 40;
  int last_decile = -1;

  // --- MCMC ---
  for (int k = 0; k < n_mcmc; ++k) {
    if (print_progress > 0) {
      double perc = 100.0 * (k + 1.0) / n_mcmc;
      int decile = (static_cast<int>(std::floor(perc)) / 10) * 10;  // 0,10,...,100
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

    // Proposal
    vec z = randn<vec>(p);
    vec beta_prop = beta + std::sqrt(ct * scale0) * (L_prop_base * z);

    double logp_prop = log_post_SL_stable(beta_prop, b0, B_inv, y, X, w, tau, L_XtWX, lambda);

    double log_alpha = std::min(0.0, logp_prop - logp_curr);
    double u = R::runif(0.0, 1.0);
    const double acc_prob = std::exp(log_alpha);

    if (u < acc_prob) {
      beta = beta_prop;
      logp_curr = logp_prop;
      ++accept;
    }

    // Scale adaptation (Robbins–Monro), target ~0.234
    double step = std::pow(k + 1.0, -0.8);
    ct = std::exp(std::log(ct) + step * (acc_prob - 0.234));
    ct = std::max(1e-4, std::min(ct, 1e4)); // clamp

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
    _["call"]        = "MCMC_BWQR_SL_STABLE"
  );
}
