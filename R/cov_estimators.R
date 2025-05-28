#' Shrinkage Estimator for Covariance Matrix
#'
#' Estimates a covariance matrix using a shrinkage approach. The sample
#' covariance matrix is shrunk toward a target matrix (the diagonal matrix of
#' sample variances) using a shrinkage intensity lambda. If lambda is not
#' provided, it is computed from the variance of the sample correlations.
#'
#' @param resid A numeric matrix of residuals.
#' @param lambda Numeric; the shrinkage intensity in [0, 1]. If missing, it is computed.
#' @param zero_mean Logical; if TRUE, residuals are assumed to have zero mean. Default is TRUE.
#'
#' @return A list containing:
#'   \item{lambda}{the optimal shrinkage intensity}
#'   \item{cov.shr}{the shrinkage covariance matrix}
#'
#' @export
shrinkage_est <- function(resid, lambda = NULL, zero_mean = TRUE){

  # Check input
  if (any(is.na(resid))) {
    print("NA presents in Residuals, function will omits rows with NA")
  }

  covm <- compute_cov_matrix(resid, zero_mean)
  target_matrix <- diag(diag(covm))

  if (is.null(lambda)) {
    r <- cov2cor(covm)
    corapn <- cov2cor(target_matrix)

    # variance of sample correlations
    v <- compute_var_cor_matrix(resid, covm, zero_mean)
    diag(v) <- 0  # diagonal variances are set to zero

    # Calculate the denominator: sum of squared differences between r and identity matrix.
    numer <- sum(v)
    denom <- sum((r - corapn)^2)

    # Compute lambda for the shrinkage estimator
    lambda <- if (denom == 0) 0 else numer / denom
  } else {
    if (lambda < 0 || lambda > 1) {
      stop("Lambda must be between 0 and 1.")
    }
  }
  lambda <- max(min(lambda, 1), 0)  # ensure lambda in [0,1]

  W <- lambda * target_matrix + (1 - lambda) * covm
  return(list(
    lambda = lambda,
    cov = W
  ))
}


#' NOVELIST Estimator for Covariance Matrix
#'
#' Estimates a covariance matrix using the NOVELIST approach. The sample
#' correlation matrix is shrunk toward its soft-thresholded version, using
#' a specified threshold (delta) and a shrinkage intensity lambda.
#' If lambda is missing, it is computed automatically from the data.
#'
#' @param resid A numeric matrix of residuals.
#' @param delta Numeric; threshold value [0, 1] applied to off-diagonal elements.
#' @param lambda Numeric; shrinkage intensity in [0, 1]. If missing, it is computed.
#' @param zero_mean Logical; if TRUE, residuals are assumed to have zero mean. Default is TRUE.
#' @param ensure_PD Logical; if TRUE, ensures the covariance matrix is positive definite.
#'
#' @return A list containing:
#'   \item{lambda}{the optimal shrinkage intensity}
#'   \item{cov.novelist}{the NOVELIST covariance matrix}
#'
#' @export
novelist_est <- function(
    resid,
    delta,
    lambda = NULL,
    zero_mean = TRUE,
    ensure_PD = TRUE
){
  # Check input
  if (any(is.na(resid))) {
    print("NA presents in Residuals, function will omits rows with NA")
  }
  if (delta < 0 || delta > 1) {
    stop("Delta must be between 0 and 1.")
  }

  T <- nrow(resid)
  p <- ncol(resid)
  covm <- compute_cov_matrix(resid, zero_mean)
  r <- cov2cor(covm)

  # Apply soft threshold to off-diagonal elements.
  r_thresh <- sign(r) * pmax(abs(r) - delta, 0)
  diag(r_thresh) <- 1 # set diagonal to 1

  # variance of sample correlations
  v <- compute_var_cor_matrix(resid, covm, zero_mean)
  diag(v) <- 0  # diagonal variances are set to zero

  # Compute numerator & denominator to compute lambda
  numer <- sum(v[abs(r) <= delta])
  denom <- sum((r - r_thresh)^2)

  if (is.null(lambda)) {
    # Compute lambda for NOVELIST using only off-diagonals (v's diagonal is zero already).
    lambda <- if (denom == 0) 0 else numer / denom
  } else {
    if (lambda < 0 || lambda > 1) {
      stop("Lambda must be between 0 and 1.")
    }
  }
  lambda <- max(min(lambda, 1), 0)  # ensure lambda in [0,1]

  # Form the NOVELIST shrunk correlation matrix.
  R_novelist <- lambda * r_thresh + (1 - lambda) * r

  # Reconstruct the covariance matrix using the NOVELIST shrunk correlation matrix.
  D_half <- diag(sqrt(diag(covm)))  # diagonal matrix of standard deviations
  W <- D_half %*% R_novelist %*% D_half

  # Ensure positive definiteness
  if (ensure_PD) {
    # Check if W is positive definite
    if (any(eigen(W, only.values = TRUE)$values <= 1e-8)) {
      W <- Matrix::nearPD(W)$mat
      # TODO: any better way?
    }
  }

  # cat("Lambda: ", lambda, "\n")
  return(list(
    lambda = lambda,
    cov = W
  ))
}


#' Compute Covariance Matrix
#'
#' Computes the covariance matrix. When residuals are assumed to have
#' zero mean, use cross-product divided by T; otherwise, cov().
#'
#' @param resid A numeric matrix of residuals with row obs and col variables.
#' @param zero_mean Logical; if TRUE, residuals are assumed to have zero mean.
#'
compute_cov_matrix <- function(resid, zero_mean) {
  T <- nrow(resid)
  if (zero_mean) {
    crossprod(stats::na.omit(resid)) / T
  } else {
    cov(stats::na.omit(resid))
  }
}


#' Compute Variance of Sample Correlation
#'
#' Computes the variance of sample correlation from residuals, which will be
#' standardised by given covariance matrix. The variance of sample correlation
#' formula differs depending on zero-mean assumption on residuals.
#'
#' @param resid A numeric matrix of residuals.
#' @param covm A covariance matrix, computed by compute_cov_matrix().
#' @param zero_mean Logical; if TRUE, residuals are assumed to have zero mean.
#'
compute_var_cor_matrix <- function(resid, covm, zero_mean) {
  T <- nrow(resid)
  # Standardise residuals
  xs <- scale(resid, center = !zero_mean, scale = sqrt(diag(covm)))
  xs <- xs[complete.cases(xs), ]  # Remove any incomplete rows
  # Estimate the variance of the sample correlations
  if (zero_mean) {
    v <- (1/(T * (T - 1))) * (crossprod(xs^2) - (1/T) * (crossprod(xs))^2)
  } else {
    v <- T/(T-1)^3 * (crossprod(xs^2) - (1/T) * (crossprod(xs))^2)
  }
}


