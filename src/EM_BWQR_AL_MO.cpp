// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::interfaces(r, cpp)]]

#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

inline double K_ratio(double z) {
  if (z < 1e-8) {
    const double z2 = z * z;
    return 1.0 + z2 / 8.0 + z2 * z2 / 192.0;
  }
  double num = R::bessel_k(z, 1.5, 1);
  double den = R::bessel_k(z, 0.5, 1);
  return num / den;
}

Rcpp::List _bwqr_weighted_em_cpp_sep(
    const Eigen::MatrixXd& y,
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& w,
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXd& Gamma,
    double tau,
    const Eigen::VectorXd& mu0,
    const Eigen::MatrixXd& sigma0,
    double a0,
    double b0,
    double eps      = 1e-6,
    int    max_iter = 1000,
    bool   verbose  = false,
    Rcpp::Nullable<Rcpp::NumericVector> fix_sigma = R_NilValue
) {

  const int n = y.rows();
  const int d = y.cols();
  const int p = x.cols();

  if (w.size() != n) stop("w must have length n.");
  if (U.rows() != d) stop("U must have d rows.");
  if (Gamma.size() > 0 && Gamma.rows() != d) stop("Gamma must have d rows (or 0 cols).");

  const int K = U.cols();
  if (K < 1) stop("U must have at least one direction (K >= 1).");

  const int G = Gamma.cols();
  if (G > 0 && (G % K != 0)) stop("Gamma.cols() must be a multiple of U.cols().");
  const int r = (K > 0) ? (G / std::max(1, K)) : 0;
  const int m_blk = p + r;

  bool prior_by_block = (mu0.size() == m_blk) && (sigma0.rows() == m_blk && sigma0.cols() == m_blk);

  MatrixXd YU = y * U;              // n x K
  MatrixXd YG; if (r > 0) YG = y * Gamma;  // n x (K*r)

  const double delta2 = 2.0 / (tau * (1.0 - tau));
  const double theta  = (1.0 - 2.0 * tau) / (tau * (1.0 - tau));

  bool use_fixed_sigma = false;
  double sigma_fixed_val = 1.0;
  if (fix_sigma.isNotNull()) {
    NumericVector fs(fix_sigma);
    if (fs.size() != 1 || !R_finite(fs[0]) || fs[0] <= 0.0)
      stop("'fix_sigma' must be a positive scalar if provided.");
    use_fixed_sigma = true;
    sigma_fixed_val = fs[0];
  }

  MatrixXd beta_prev = MatrixXd::Zero(K, m_blk);  // K x (p+r)
  VectorXd sigma_prev = VectorXd::Constant(K, 1.0);

  for (int k = 0; k < K; ++k) {
    MatrixXd Xstar(n, m_blk);
    Xstar.leftCols(p) = x;
    if (r > 0) {
      Xstar.rightCols(r) = YG.block(0, k * r, n, r);
    }
    VectorXd yk = YU.col(k);
    // OLS
    VectorXd bk = (Xstar.transpose() * Xstar).ldlt().solve(Xstar.transpose() * yk);
    beta_prev.row(k) = bk.transpose();
    if (use_fixed_sigma) {
      sigma_prev(k) = sigma_fixed_val;
    } else {
      VectorXd res = yk - Xstar * bk;
      double s2 = (res.array().square() * w.array()).mean();
      sigma_prev(k) = std::max(1e-3, s2);
    }
  }

  // ===== EM =====
  MatrixXd beta_curr = beta_prev;
  VectorXd sigma_curr = sigma_prev;

  for (int iter = 0; iter < max_iter; ++iter) {

    for (int k = 0; k < K; ++k) {
      const double inv_ds = 1.0 / (delta2 * sigma_prev(k));

      MatrixXd Xstar(n, m_blk);
      Xstar.leftCols(p) = x;
      if (r > 0) {
        Xstar.rightCols(r) = YG.block(0, k * r, n, r);
      }

      VectorXd y_aux_k(n), e_nu_inv_k(n), e_nu_k(n);
      for (int i = 0; i < n; ++i) {
        double fit = Xstar.row(i).dot(beta_prev.row(k));
        const double rres = YU(i, k) - fit;

        const double as = std::max(w(i) * rres * rres * inv_ds, 1e-12);
        const double bs = std::max(w(i) * (2.0 / sigma_prev(k) + theta * theta * inv_ds), 1e-12);

        const double sa = std::sqrt(as);
        const double sb = std::sqrt(bs);

        e_nu_inv_k(i) = sb / sa;
        const double aux_obj = sa * sb;
        e_nu_k(i)     = (1.0 / e_nu_inv_k(i)) * K_ratio(aux_obj);
        y_aux_k(i)    = YU(i, k) - theta / e_nu_inv_k(i);
      }
      if (!e_nu_inv_k.allFinite() || !e_nu_k.allFinite())
        stop("NaN/Inf in E-step (direction k) - check data or algorithm.");

      VectorXd mu0_blk(m_blk);
      MatrixXd S0inv_blk(m_blk, m_blk);
      if (prior_by_block) {
        mu0_blk = mu0;
        S0inv_blk = sigma0.ldlt().solve(MatrixXd::Identity(m_blk, m_blk));
      } else {
        mu0_blk = mu0.segment(k * m_blk, m_blk);
        S0inv_blk = sigma0.block(k * m_blk, k * m_blk, m_blk, m_blk)
                          .ldlt().solve(MatrixXd::Identity(m_blk, m_blk));
      }

      MatrixXd A = S0inv_blk;
      VectorXd B = S0inv_blk * mu0_blk;

      for (int i = 0; i < n; ++i) {
        const double coef = w(i) * e_nu_inv_k(i) * inv_ds;
        A.noalias() += coef * (Xstar.row(i).transpose() * Xstar.row(i));
        B.noalias() += coef * Xstar.row(i).transpose() * y_aux_k(i);
      }

      A.diagonal().array() += 1e-8;
      VectorXd beta_k = A.ldlt().solve(B);
      beta_curr.row(k) = beta_k.transpose();

      if (use_fixed_sigma) {
        sigma_curr(k) = sigma_fixed_val;
      } else {
        double term1 = 0.0, term2 = 0.0, term3 = 0.0;
        for (int i = 0; i < n; ++i) {
          const double r2 = y_aux_k(i) - Xstar.row(i).dot(beta_k);
          term1 += w(i) * e_nu_inv_k(i) * r2 * r2;
          term2 += w(i) * theta * theta * (e_nu_k(i) - 1.0 / e_nu_inv_k(i));
          term3 += w(i) * e_nu_k(i);
        }
        term1 /= (2.0 * delta2);
        term2 /= (2.0 * delta2);
        const double denom = (3.0 * n + a0 + 1.0) / 2.0;

        sigma_curr(k) = (term1 + term2 + term3 + b0) / denom;
        if (!std::isfinite(sigma_curr(k)))
          stop("sigma_k became NaN/Inf - algorithm diverged.");
      }
    }

    double diff = 0.0;
    diff += (beta_curr - beta_prev).cwiseAbs().sum();
    diff += (sigma_curr - sigma_prev).cwiseAbs().sum();

    if (verbose && (iter % 50 == 0)) {
      Rcpp::Rcout << "iter=" << iter << "  diff=" << diff
                  << "  mean(sigma)=" << sigma_curr.mean() << std::endl;
    }

    if (diff < eps) {
      return List::create(
        _["beta"]      = beta_curr,   // K x (p+r)
        _["sigma"]     = sigma_curr,  // K
        _["iter"]      = iter + 1,
        _["converged"] = true
      );
    }

    beta_prev  = beta_curr;
    sigma_prev = sigma_curr;
  }

  if (verbose)
    Rcpp::Rcout << "EM reached max_iter without convergence (diff > eps)." << std::endl;

  return List::create(
    _["beta"]      = beta_prev,
    _["sigma"]     = sigma_prev,
    _["iter"]      = max_iter,
    _["converged"] = false
  );
}
