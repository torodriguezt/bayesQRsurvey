// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <R_ext/Utils.h>
#include <R_ext/Print.h>

using namespace Rcpp;
using namespace arma;

constexpr double PI_CONST = 3.14159265358979323846;

inline bool safe_chol(arma::mat& L, const arma::mat& A, double& ridge,
                      const bool lower = true) {
  const double diag_mean = arma::mean(A.diag());
  ridge = std::max(ridge, 1e-12 * std::max(1.0, std::abs(diag_mean)));
  arma::mat M;
  for (int i = 0; i < 10; ++i) {
    M = A + ridge * arma::eye<arma::mat>(A.n_rows, A.n_cols);
    if (arma::chol(L, M, lower ? "lower" : "upper")) return true;
    ridge *= 10.0;
  }
  return false;
}

inline arma::vec spd_solve_from_chol(const arma::mat& L, const arma::vec& s) {
  arma::vec v = arma::solve(arma::trimatl(L), s, arma::solve_opts::fast);
  arma::vec z = arma::solve(arma::trimatu(L.t()), v, arma::solve_opts::fast);
  return z;
}

static double log_post(const arma::vec& beta,
                       const arma::vec& b0,
                       const arma::mat& B_inv,
                       const arma::vec& y,
                       const arma::mat& X,
                       const arma::vec& w,
                       const arma::vec& d0,
                       const arma::mat& Cpair,
                       bool have_pi2,
                       double tau)
{
  const arma::uword p = beta.n_elem;

  arma::vec diff = beta - b0;
  double lp = -0.5 * arma::dot(diff, B_inv * diff);

  arma::vec res = y - X * beta;
  arma::vec psi = tau - arma::conv_to<arma::vec>::from(res < 0.0);

  arma::vec U_S = X.t() * (w % psi);

  arma::mat Omega(p, p);
  if (have_pi2) {
    arma::mat Xtil = X.each_col() % psi;
    Omega = Xtil.t() * (Cpair * Xtil);
  } else {
    arma::vec s = arma::sqrt(d0) % arma::abs(psi);
    arma::mat Xs = X.each_col() % s;
    Omega = Xs.t() * Xs;
  }

  arma::mat L;
  double ridge = 0.0;
  if (!safe_chol(L, Omega, ridge, true))
    return -std::numeric_limits<double>::infinity();

  arma::vec z   = spd_solve_from_chol(L, U_S);
  double quad   = arma::dot(U_S, z);
  double logdet = 2.0 * arma::accu(arma::log(L.diag()));

  lp += -0.5 * quad
        - 0.5 * (static_cast<double>(p) * std::log(2.0 * PI_CONST) + logdet);

  return lp;
}

Rcpp::List _mcmc_bwqr_ap_cpp(const arma::vec& y,
                             const arma::mat& X,
                             const arma::vec& w,
                             int n_mcmc,
                             int burnin,
                             int thin,
                             double tau = 0.5,
                             Rcpp::Nullable<Rcpp::NumericVector> b_prior_mean = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericMatrix> B_prior_prec = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericMatrix> pi_matrix = R_NilValue,
                             int print_progress = 0)
{
  if (y.n_elem != X.n_rows || w.n_elem != y.n_elem)
    stop("Dimensions of y, X and w must match.");
  if (!(tau > 0.0 && tau < 1.0))
    stop("tau must be in (0,1).");
  if (n_mcmc <= 0)  stop("n_mcmc must be positive");
  if (burnin < 0)   stop("burnin must be non-negative");
  if (thin   <= 0)  stop("thin must be positive");
  if (burnin >= n_mcmc)
    stop("burnin must be less than n_mcmc");
  if (w.min() <= 0.0) stop("All weights must be > 0.");

  const int p = X.n_cols;
  const int n = y.n_elem;

  arma::vec b0 = b_prior_mean.isNotNull() ? Rcpp::as<arma::vec>(b_prior_mean) : arma::zeros<vec>(p);
  if ((int)b0.n_elem != p)
    stop("b_prior_mean must have length equal to ncol(X)");

  arma::mat B_inv;
  if (B_prior_prec.isNotNull()) {
    B_inv = Rcpp::as<arma::mat>(B_prior_prec);
    if ((int)B_inv.n_rows != p || (int)B_inv.n_cols != p)
      stop("B_prior_prec must be a pxp matrix");
  } else {
    B_inv = arma::eye<mat>(p, p) / 100.0;
  }

  arma::vec pi1(n);
  arma::vec d0(n);
  arma::mat Cpair;
  bool have_pi2 = false;

  if (pi_matrix.isNotNull()) {
    arma::mat Pm = Rcpp::as<arma::mat>(pi_matrix);
    if ((int)Pm.n_rows != n || (int)Pm.n_cols != n)
      stop("pi_matrix must be an n x n matrix (n = length(y)).");

    pi1 = Pm.diag();
    if (pi1.min() <= 0.0 || pi1.max() > 1.0)
      stop("First-order inclusion probabilities (diag(pi_matrix)) must be in (0,1].");

    arma::mat OffAbs = arma::abs(Pm);
    OffAbs.diag().zeros();
    have_pi2 = (OffAbs.max() > 0.0);

    if (have_pi2) {
      Cpair.zeros(n, n);
      Cpair.diag() = (1.0 - pi1) / arma::square(pi1);
      for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
          double pij = Pm(i, j);
          if (pij <= 0.0) pij = Pm(j, i);
          if (pij <= 0.0)
            stop("Second-order inclusion probabilities (pi_ij) must be > 0.");
          double c = (pij - pi1(i) * pi1(j)) / (pij * pi1(i) * pi1(j));
          Cpair(i, j) = c;
          Cpair(j, i) = c;
        }
      }
    }
  } else {
    pi1 = 1.0 / w;
    if (pi1.max() > 1.0)
      Rcpp::warning("Some implied inclusion probabilities 1/w_i exceed 1; check the weights.");
  }

  d0 = (1.0 - pi1) / arma::square(pi1);

  if (!have_pi2) {
    Rcpp::warning(
      "Second-order inclusion probabilities not provided: using a biased "
      "Horvitz-Thompson variance estimator (first term only). Supply "
      "'pi_matrix' (diagonal = first-order, triangle = second-order) for the "
      "unbiased estimator.");
  }

  arma::vec w2 = w % w;
  arma::mat A_prop = (1.0 / n) * (X.t() * (X.each_col() % w2));
  arma::mat L_A;
  double ridge_A = 0.0;
  if (!safe_chol(L_A, A_prop, ridge_A, true))
    stop("Cholesky of (1/n) sum w_i^2 x_i x_i' failed even with ridge.");
  arma::mat I_p = arma::eye<mat>(p, p);
  arma::mat L_prop_base = arma::solve(arma::trimatu(L_A.t()), I_p, arma::solve_opts::fast);

  const int n_keep = (n_mcmc - burnin) / thin;
  arma::mat beta_out((n_keep > 0 ? n_keep : 0), p, fill::none);
  int accept = 0, k_out = 0;

  arma::vec beta = arma::solve(X, y);
  double ct = 1.0;

  double logp_curr = log_post(beta, b0, B_inv, y, X, w, d0, Cpair, have_pi2, tau);

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

    arma::vec zr = arma::randn<vec>(p);
    arma::vec beta_prop = beta + std::sqrt(ct * tau * (1.0 - tau)) * (L_prop_base * zr);

    double logp_prop = log_post(beta_prop, b0, B_inv, y, X, w, d0, Cpair, have_pi2, tau);

    double log_acc  = std::min(0.0, logp_prop - logp_curr);
    double acc_prob = std::exp(log_acc);

    if (R::runif(0.0, 1.0) < acc_prob) {
      beta = beta_prop;
      logp_curr = logp_prop;
      ++accept;
    }

    double step = std::pow(k + 1.0, -0.8);
    ct = std::exp(std::log(ct) + step * (acc_prob - 0.234));
    ct = std::max(1e-4, std::min(ct, 1e4));

    if (k >= burnin && ((k - burnin) % thin == 0)) {
      if (k_out < n_keep) beta_out.row(k_out++) = beta.t();
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
    _["n_mcmc"]      = n_mcmc,
    _["burnin"]      = burnin,
    _["thin"]        = thin,
    _["n_samples"]   = n_keep,
    _["call"]        = "MCMC_BWQR_AP"
  );
}
