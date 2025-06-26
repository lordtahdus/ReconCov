// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Load novelist_est_cpp from novelist.cpp
Rcpp::List novelist_est_cpp(const arma::mat& resid,
                            double delta,
                            Rcpp::Nullable<double> lambda_in = R_NilValue,
                            bool zero_mean = true);

//' Enforce positive definiteness on a matrix in C++
//'
//' This function performing an eigen decomposition on a symmetric matrix,
//' thresholding the eigenvalues, and reconstructing the matrix.
//'
//' @param W The input matrix to be made positive definite.
//' @param threshold A small value to replace negative eigenvalues (1e-8 default).
//' @return A positive definite matrix.
//'
//' @export
// [[Rcpp::export]]
arma::mat make_PD_cpp(
    const arma::mat& W,
    double tol = 1e-8
) {
  arma::vec eigval;
  arma::mat eigvec;

  if (!arma::eig_sym(eigval, eigvec, W)) {
    Rcpp::stop("Eigen decomposition failed.");
  }
  eigval.transform([&](double x) {
    return (x < tol) ? tol : x;
  });
  arma::mat W_PD = eigvec * diagmat(eigval) * eigvec.t();

  // Symmetrize to fix numeric asymmetry
  W_PD = 0.5 * (W_PD + W_PD.t());
  return W_PD;
}


//' Reconcile forecasts using MinT formula in C++
//'
//' @seealso \code{\link[=reconcile_mint]}
//' @export
// [[Rcpp::export]]
arma::mat reconcile_mint_cpp(
    const arma::mat& base_forecasts,
    const arma::mat& S,
    const arma::mat& W
) {
  arma::mat R = S.t() *  solve(W);
  arma::mat P = solve(R * S) * R;

  // if (base_forecasts.n_rows == 1 || base_forecasts.n_cols == 1) {
  //   arma::vec fc_vec = vectorise(base_forecasts);
  //   arma::mat out = trans(S * P * fc_vec);
  //   return out;
  // } else if (S.n_rows == base_forecasts.n_cols) {
  //   return trans(S * P * trans(base_forecasts));
  // } else {
  //   Rcpp::stop("Dimensions do not conform: base_forecasts must be vector or T x p matrix.");
  // }
  return trans(S * P * base_forecasts.t()).t();
}


//' Outsource the two for loops in novelist_cv into C++
//'
// [[Rcpp::export]]
Rcpp::List novelist_cov_grid_cpp(
    const arma::mat& resid,
    const arma::vec& deltas,
    int window_size,
    bool zero_mean = true
) {
  // basic checks
  if (window_size < 1 || window_size >= static_cast<int>(resid.n_rows))
    Rcpp::stop("window_size out of range");
  if (deltas.n_elem == 0)
    Rcpp::stop("deltas must be non‑empty");

  int T      = resid.n_rows;
  int p      = resid.n_cols;   // not used but sanity
  int n_win  = T - window_size;
  int n_del  = deltas.n_elem;

  Rcpp::List win_list(n_win);  // initialise list of windows

  for (int k = 0; k < n_win; ++k) {
    // rows [k, k+window_size-1]
    arma::mat train = resid.rows(k, k + window_size - 1);
    Rcpp::List cov_list(n_del);    // initialise list of covariances for this window

    for (int j = 0; j < n_del; ++j) {
      double delta = deltas[j];
      arma::mat covm = Rcpp::as<arma::mat>(
        novelist_est_cpp(train, delta, R_NilValue, zero_mean)["cov"]
      );
      cov_list[j] = covm;  // store covariance matrix
    }
    win_list[k] = cov_list;
  }

  win_list.attr("deltas") = deltas;
  return win_list;
}

//' Rolling Cross-Validation for NOVELIST Threshold Selection (C++ accelerated)
//'
//' @seealso \code{\link[=novelist_cv]}
//' @note This function use different enforcing PD method.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List novelist_cv_cpp(
    const arma::mat&  y,
    const arma::mat&  y_hat,
    const arma::mat&  S,
    int               window_size,
    const arma::vec&  deltas,
    bool              zero_mean = true,
    bool              ensure_PD = true,
    double            PD_tol    = 1e-8
) {
  // Checks
  if (y.n_rows != y_hat.n_rows || y.n_cols != y_hat.n_cols)
    stop("Dimensions of y and y_hat must match.");
  if (window_size <= 1 || window_size >= static_cast<int>(y.n_rows))
    stop("Invalid window_size.");
  if (deltas.n_elem == 0)
    stop("deltas must be non-empty.");

  int T     = y.n_rows;
  int n_win = T - window_size;
  int n_del = deltas.n_elem;

  arma::mat resid = y - y_hat;
  arma::mat cv_err(n_win, n_del, fill::zeros);

  // Rolling CV
  for (int k = 0; k < n_win; ++k) {

    arma::mat train = resid.rows(k, k + window_size - 1);

    int v_idx = k + window_size;               // validation index
    arma::rowvec actual_next = y.row(v_idx);
    arma::rowvec fitted_next = y_hat.row(v_idx);

    for (int j = 0; j < n_del; ++j) {
      // novelist
      arma::mat W = Rcpp::as<arma::mat>(
        novelist_est_cpp(train, deltas[j], R_NilValue, zero_mean)["cov"]
      );
      if (ensure_PD) W = make_PD_cpp(W, PD_tol);

      // Reconcile & MSE
      arma::rowvec recon = reconcile_mint_cpp(fitted_next, S, W);
      cv_err(k, j)      = mean(square(actual_next - recon));
    }
  }

  // Aggregation & best δ
  arma::rowvec mean_err = mean(cv_err, 0);
  uword best_j; mean_err.min(best_j);
  double delta_star = deltas(best_j);

  // Refit NOVELIST on *all* residuals with best δ
  List final_est = novelist_est_cpp(resid, delta_star, R_NilValue, zero_mean);
  arma::mat cov_final = final_est["cov"];
  double lambda_star  = final_est["lambda"];
  if (ensure_PD) cov_final = make_PD_cpp(cov_final, PD_tol);

  // Return
  return List::create(
    _["delta"]           = delta_star,
    _["lambda"]          = lambda_star,
    _["cov"]             = cov_final,
    _["errors_by_delta"] = Rcpp::NumericVector(mean_err.begin(), mean_err.end())
  );
}
