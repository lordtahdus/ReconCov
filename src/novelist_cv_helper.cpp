// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


Rcpp::List novelist_est_cpp(const arma::mat& resid,
                            double delta,
                            Rcpp::Nullable<double> lambda_in = R_NilValue,
                            bool zero_mean = true);


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

