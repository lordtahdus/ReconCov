// novelist.cpp — RcppArmadillo implementation of the NOVELIST covariance estimator
//
// Compile in R with something like:
//   Rcpp::sourceCpp("novelist.cpp")
// then call `novelist_est_cpp(resid, delta, R_NilValue, true, true)`
//
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//---------------------------------------------------------------
// Helpers
//---------------------------------------------------------------

// fast covariance computation (zero‑mean or sample mean)
arma::mat compute_cov_matrix_cpp(
    const arma::mat& resid,
    bool zero_mean
) {
  const double T = static_cast<double>(resid.n_rows);

  if (zero_mean) {
    return (resid.t() * resid) / T;
  }
  return arma::cov(resid, 0); // use sample covariance
}

// variance of sample correlations (Schäfer–Strimmer 2005)
arma::mat compute_var_cor_matrix_cpp(
    const arma::mat& resid,
    const arma::mat& covm,
    bool zero_mean
) {
  const double T = static_cast<double>(resid.n_rows);
  const arma::rowvec sd = sqrt(covm.diag()).t();

  // Standardise residuals (optionally centre first)
  arma::mat xs = resid;
  if (!zero_mean) {
    // De-mean
    arma::rowvec mu = mean(xs, 0);
    xs.each_row() -= mu;
  }
  xs.each_row() /= sd;

  // xs^2
  arma::mat xs2 = xs % xs;

  arma::mat cp_xs2 = xs2.t() * xs2;        // crossprod(xs^2)
  arma::mat cp_xs = xs.t() * xs;           // crossprod(xs)
  arma::mat cp_xs_sq = cp_xs % cp_xs;      // (crossprod(xs))^2

  arma::mat v;
  if (zero_mean) {
    v = (1.0 / (T * (T - 1.0))) * (cp_xs2 - (1.0 / T) * cp_xs_sq);
  } else {
    v = (T / std::pow(T - 1.0, 3)) * (cp_xs2 - (1.0 / T) * cp_xs_sq);
  }
  return v;
}


//---------------------------------------------------------------
// Main NOVELIST estimator
//---------------------------------------------------------------

//' @title NOVELIST covariance estimator, Rcpp version
//'
//' @seealso \code{\link[=novelist_est]}
//' @useDynLib ReconCov
//' @import Rcpp
//' @importFrom Rcpp sourceCpp
//' @import RcppArmadillo
//' @export
// [[Rcpp::export]]
Rcpp::List novelist_est_cpp(
    const arma::mat& resid,
    double delta,
    Rcpp::Nullable<double> lambda_in = R_NilValue,
    bool zero_mean  = true
) {

  if (delta < 0.0 || delta > 1.0)
    stop("delta must be between 0 and 1");

  // Sample covariance & correlation
  arma::mat covm = compute_cov_matrix_cpp(resid, zero_mean);

  arma::vec sd = sqrt(covm.diag());
  arma::mat r = covm;
  r.each_col() /= sd;     // divide by sd_i (columns)
  r.each_row() /= sd.t(); // divide by sd_j (rows)

  // Soft‑threshold correlation matrix
  arma::mat r_thresh = r;
  for (uword j = 0; j < r.n_cols; ++j) {      // Using uword (rather than plain int) avoids any signed/unsigned mismatch when you compare against r.n_cols (which itself is a uword).
    for (uword i = 0; i < r.n_rows; ++i) {
      if (i == j) {
        r_thresh(i, j) = 1.0;           // keep diag at 1
      } else {
        double val = r(i, j);
        double thr_val = std::max(std::abs(val) - delta, 0.0);
        r_thresh(i, j) = (val >= 0.0 ? thr_val : -thr_val);   // sign
      }
    }
  }

  // Variance of correlations & lambda
  arma::mat v = compute_var_cor_matrix_cpp(resid, covm, zero_mean);
  v.diag().zeros();             // diag variance = 0

  double numer = 0.0, denom = 0.0;
  for (uword j = 1; j < r.n_cols; ++j) {
    for (uword i = 0; i < j; ++i) {          // off‑diagonals only
      if (std::abs(r(i, j)) <= delta)
        numer += v(i, j) + v(j, i);      // symmetric sum
      double diff = r(i, j) - r_thresh(i, j);
      denom += diff * diff * 2.0;          // both triangles
    }
  }

  double lambda;
  if (lambda_in.isNotNull()) {
    lambda = Rcpp::as<double>(lambda_in);
    if (lambda < 0.0 || lambda > 1.0)
      stop("lambda must be in [0,1]");
  } else {
    lambda = denom == 0.0 ? 0.0 : numer / denom; // avoid /0
  }
  lambda = std::max(0.0, std::min(1.0, lambda));

  // Shrunk correlation & covariance
  arma::mat R_novelist = lambda * r_thresh + (1.0 - lambda) * r;
  arma::mat D_half = diagmat(sd);
  arma::mat W = D_half * R_novelist * D_half;

  // TODO: ensure positive definiteness

  // Return list equivalent to R version
  return Rcpp::List::create(
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("cov")   = W
  );
}
