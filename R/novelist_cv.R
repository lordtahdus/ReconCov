#' Rolling Cross-Validation for NOVELIST Threshold Selection
#'
#' This function performs a rolling window approach to select the threshold 'delta'
#' for the NOVELIST covariance estimator by minimizing some measure of forecast error
#' on validation windows.
#'
#' @param y A matrix of actual values (dim T by p),
#' @param y_hat A matrix of fitted values (dim T by p).
#' @param S A matrix of reconciliation structure (dim p by b; b is number of bottom series).
#' @param window_size The length of the rolling window n.
#' @param deltas A numeric vector of candidate threshold values in [0,1].
#' @param reconcile_forecast .....
#' @param zero_mean Logical, whether to treat the residuals as zero mean in the covariance.
#' @param error_metric An error measure function for given actual and reconciled forecasts.
#' @param ensure_PD Logical, whether to ensure the covariance matrix is positive definite.
#'
#' @return A list containing:
#'   \item{delta}{optimal threshold}
#'   \item{lambda}{optimal shrinkage intensity}
#'   \item{cov.novelist}{optimal NOVELIST covariance matrix}
#'   \item{errors_by_delta}{vector of average validation error for each delta}
#'
#' @details
#' This approach loops through each possible validation segment in a rolling manner:
#'   for i in window_size:(T-1), use resid[(i-window_size+1):i, ] as training data
#'   then produce a reconciled forecast for time i+1 and measure the error.
#' The final best delta is the one that yields the lowest error measure.
#'
#' @export
novelist_cv <- function(
    y,
    y_hat,
    S,
    window_size,
    deltas = seq(0, 1, by = 0.05),
    # reconcile_forecast = NULL,
    zero_mean = TRUE,
    error_metric = function(actual, fc) mean((actual - fc)^2, na.rm = TRUE),
    ensure_PD = TRUE,
    message = TRUE
) {
  T <- nrow(y)
  p <- ncol(y)

  # Preliminary checks
  if(any(dim(y) != dim(y_hat))) {
    stop("Dimensions of 'y' and 'y_hat' must match.")
  }
  if(T <= window_size || window_size <= 1 || (window_size %% 1 != 0)) {
    stop("Window size must be integer, greater than 1 and less than nrow.")
  }

  # Calculate base residuals
  resid <- y - y_hat
  # Store errors
  cv_errors <- matrix(NA, nrow = T - window_size, ncol = length(deltas))

  # Rolling over each possible validation step
  # i means the "last index" of the training set is i.
  # Then the "validation" point is i+1.
  for (i in window_size:(T - 1)) {

    # Training residuals from (i-window_size+1) to i
    train_resid <- resid[(i - window_size + 1):i, , drop=FALSE] # drop=F keep matrix structure

    # The actual data at time i+1
    actual_next <- y[i+1, ]
    fitted_next <- y_hat[i+1, ]

    # Now loop over candidate threshold delta
    for(delta in deltas){

      # Estimate a NOVELIST covariance using the train_resid
      cov_novelist <- novelist_est(
        resid     = train_resid,
        delta     = delta,
        zero_mean = zero_mean,
        ensure_PD = ensure_PD
      )$cov

      if (any(eigen(cov_novelist)$values <= 1e-12)) {
        # stop(sort(eigen(cov_novelist)$values)[1])
        stop("The covariance matrix is not positive definite, cannot reconcile. Try ensure_PD = T")
      }

      recon_fc <- fitted_next
      recon_fc <- reconcile_mint(
        base_forecasts = fitted_next,
        S = S,
        W = cov_novelist
      )

      # Compute error measure for time i+1:
      err_val <- error_metric(actual_next, recon_fc)

      # store
      idx_i <- i - window_size + 1
      idx_d <- which(deltas == delta)
      cv_errors[idx_i, idx_d] <- err_val
    } # end loop over deltas

    # TODO: Refine Progress Printing
    if (message) {
      if (idx_i %% 5 == 0 || idx_i == (T - window_size)) {
        cat("Iteration: ", idx_i, " / ", (T - window_size), "\n")
      }
    }
  }
  cat("\n")

  # Now compute average error across all rolling iterations for each delta
  mean_errors <- colMeans(cv_errors, na.rm = TRUE)
  names(mean_errors) <- deltas

  # find best delta
  best_idx <- which.min(mean_errors)
  delta_star <- deltas[best_idx]

  # Refit final NOVELIST covariance using entire resid
  novelist_results <- novelist_est(
    resid     = resid,
    delta     = delta_star,
    zero_mean = zero_mean
  )

  final_cov_novelist <- novelist_results$cov
  lambda_star <- novelist_results$lambda

  # Check final positive definiteness
  if (any(eigen(final_cov_novelist)$values <= 1e-8)) {
    stop("The final covariance matrix is not positive definite.")
  }

  # Return results
  return(list(
    delta = delta_star,
    lambda = lambda_star,
    cov = final_cov_novelist,
    errors_by_delta = mean_errors
  ))
}


#' Rolling Cross-Validation for NOVELIST Threshold Selection (C++ grid)
#'
#' This function uses Rcpp to speed up the rolling window computation of
#' NOVELIST covariance matrices across a grid of threshold values. It avoids
#' calling nearPD or doing any matrix inversion inside C++.
#'
#' @inheritParams novelist_cv
#' @return A list containing:
#'   \item{delta}{optimal threshold}
#'   \item{lambda}{optimal shrinkage intensity}
#'   \item{cov.novelist}{optimal NOVELIST covariance matrix}
#'   \item{errors_by_delta}{vector of average validation error for each delta}
#'
#' @export
novelist_cv_grid <- function(
    y,
    y_hat,
    S,
    window_size,
    deltas = seq(0, 1, by = 0.05),
    zero_mean = TRUE,
    error_metric = function(actual, fc) mean((actual - fc)^2, na.rm = TRUE),
    ensure_PD = TRUE
) {
  T <- nrow(y)
  p <- ncol(y)

  if (any(dim(y) != dim(y_hat))) {
    stop("Dimensions of 'y' and 'y_hat' must match.")
  }
  if (T <= window_size || window_size <= 1 || (window_size %% 1 != 0)) {
    stop("Window size must be integer, greater than 1 and less than nrow.")
  }

  resid <- y - y_hat
  deltas <- as.numeric(deltas)

  # 1. C++ accelerated computation of all covariances
  grid <- novelist_cov_grid_cpp(resid, deltas, window_size, zero_mean)

  n_win <- length(grid)
  n_del <- length(deltas)

  # 2. Error matrix
  cv_errors <- matrix(NA_real_, nrow = n_win, ncol = n_del)

  # Loop only ONCE over rolling windows (each row of errors)
  for (k in seq_len(n_win)) {
    i <- window_size + k
    actual_next <- y[i, ]
    fitted_next <- y_hat[i, ]
    cov_list <- grid[[k]]

    # Vectorized over deltas using vapply
    cv_errors[k, ] <- vapply(seq_along(cov_list), function(j) {
      W <- cov_list[[j]]
      if (ensure_PD) {
        W <- as.matrix(Matrix::nearPD(W)$mat)
      }
      recon <- reconcile_mint(fitted_next, S, W)
      error_metric(actual_next, recon)
    }, numeric(1))
  }

  mean_errors <- colMeans(cv_errors, na.rm = TRUE)
  names(mean_errors) <- deltas
  best_idx <- which.min(mean_errors)
  delta_star <- deltas[best_idx]

  # Final re-estimation
  final_result <- novelist_est(
    resid = resid,
    delta = delta_star,
    zero_mean = zero_mean,
    ensure_PD = ensure_PD
  )

  return(list(
    delta = delta_star,
    lambda = final_result$lambda,
    cov = final_result$cov,
    errors_by_delta = mean_errors
  ))
}
