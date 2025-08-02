## Covariance Estimators for Forecast Reconciliation (beta)
*Note: R Package for Research Purpose*

This package implements Min Trace reconciliation method with alternative covariance estimators.

- Shrinkage - the original estimator used in Min Trace paper
- NOVELIST - NOVEL Integration of the Sample and Thresholded Covariance estimator
- ...

For NOVELIST estimator, there is a cross validation algorithm to select the optimal thresholding value (no closed-form expression), which also has its C++ version to speed up the process (50% faster than R version).

This package is also built for simulations: generating hierarchical data with customisable VAR(1) and correlation data-generating-process structures.

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
