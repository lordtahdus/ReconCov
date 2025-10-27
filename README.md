# Covariance Estimators for Minimum Trace Forecast Reconciliation (beta)
*Note: R Package for Research Purpose*

This package implements Min Trace reconciliation method with alternative covariance estimators.

For more details, please refer to the main thesis document [here](https://github.com/lordtahdus/Recon_Honours_Thesis/blob/master/Honours_thesis_25_VincentSu.pdf).

This package is also built for simulations: hierarchical data with VAR(1) coefficients with realistic innovation process.

## Main Functions

`reconcile_mint()` : Reconcile base forecasts using Min Trace method with covariance estimate.

`shrinkage_est()` : Compute shrinkage covariance estimator.
`novelist_est()` : Compute NOVELIST covariance estimator.
`novelist_cv()` : Cross validation for NOVELIST threshold selection.

`shrinkage_pc_est()` : Compute PC-adjusted shrinkage covariance estimator.
`novelist_pc_est()` : Compute PC-adjusted NOVELIST covariance estimator.
`novelist_pc_cv()` : Cross validation for PC-adjusted NOVELIST threshold selection.

**C++ Versions:**

Below implementations speed up computation by 50%.

`reconcile_mint_cpp()` : C++ version of Min Trace reconciliation.
`novelist_est_cpp()` : C++ version of NOVELIST covariance estimator.
`novelist_cv_cpp()` : C++ version of cross validation for NOVELIST threshold selection.

All functions have detailed documentation accessible via R's help system, e.g., `?reconcile_mint` or `help(reconcile_mint)`.

## Installation

You can install the latest development release directly from GitHub:

```r
# install 'remotes' if you don’t already have it
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# install ReconCov from GitHub
remotes::install_github("your-username/ReconCov")
```

Loading the Package:

```r
library(ReconCov)
```


<br> <br> <br> <br>

For working paper markdown production, click [here](https://github.com/lordtahdus/Recon_Honours_Thesis) (private)
