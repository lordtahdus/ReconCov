#' Rolling Cross-Validation for NOVELIST Threshold Selection
#'
#' This function performs a rolling window approach to select the threshold 'delta'
#' for the NOVELIST covariance estimator by minimizing some measure of forecast error
#' on validation windows.
#'
#' @param y A matrix of actual values (dim T by p),
#' @param base_forecasts A matrix of base forecasts (dim T by p).
#' @param S A matrix of reconciliation structure (dim p by b; b is number of bottom series).
#' @param window_size The length of the rolling window n.
#' @param deltas A numeric vector of candidate threshold values in [0,1].
#' @param reconcile_forecast .....
#' @param zero_mean Logical, whether to treat the residuals as zero mean in the covariance.
#' @param error_metric An error measure function for given actual and reconciled forecasts.
#'
#' @return A list containing:
#'   \item{delta}{optimal threshold}
#'   \item{lambda}{optimal shrinkage intensity}
#'   \item{cov.novelist}{optimal NOVELIST covariance matrix}
#'   \item{errors_by_delta}{vector of average validation error for each delta}
#'
#' @details
#' This approach loops through each possible validation segment in a rolling manner:
#'   for i in window_size:(T-1), use base_resid[(i-window_size+1):i, ] as training data
#'   then produce a reconciled forecast for time i+1 and measure the error.
#' The final best delta is the one that yields the lowest error measure.
#'
#' @export
novelist_cv <- function(
    y,
    base_forecasts,
    S,
    window_size,
    deltas = seq(0, 1, by = 0.05),
    # reconcile_forecast = NULL,
    zero_mean = TRUE,
    error_metric = function(actual, fc) mean((actual - fc)^2, na.rm = TRUE)
) {
  T <- nrow(y)
  p <- ncol(y)

  # Preliminary checks
  if(any(dim(y) != dim(base_forecasts))) {
    stop("Dimensions of 'y' and 'base_forecasts' must match.")
  }
  if(T <= window_size || window_size <= 1 || (window_size %% 1 != 0)) {
    stop("Window size must be integer, greater than 1 and less than nrow.")
  }

  # Calculate base residuals
  base_resid <- y - base_forecasts
  # Store errors
  cv_errors <- matrix(NA, nrow = T - window_size, ncol = length(deltas))

  # Rolling over each possible validation step
  # i means the "last index" of the training set is i.
  # Then the "validation" point is i+1.
  for (i in window_size:(T - 1)) {

    # Training residuals from (i-window_size+1) to i
    train_resid <- base_resid[(i - window_size + 1):i, , drop=FALSE] # drop=F keep matrix structure

    # The actual data at time i+1
    actual_next <- y[i+1, ]
    base_fc_next <- base_forecasts[i+1, ]

    # Now loop over candidate threshold delta
    for(delta in deltas){

      # Estimate a NOVELIST covariance using the train_resid
      cov_novelist <- novelist_est(
        resid     = train_resid,
        delta     = delta,
        zero_mean = zero_mean
      )$cov.novelist

      if (any(eigen(cov_novelist)$values <= 1e-8)) {
        stop("The covariance matrix is not positive definite. Function will be fixed")
      }

      recon_fc <- base_fc_next
      recon_fc <- reconcile_mint(
        base_forecasts = base_fc_next,
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
    if (idx_i %% 5 == 0 || idx_i == (T - window_size)) {
      cat("Iteration: ", idx_i, " / ", (T - window_size), "\n")
    }
  }
  cat("\n")

  # Now compute average error across all rolling iterations for each delta
  mean_errors <- colMeans(cv_errors, na.rm = TRUE)
  names(mean_errors) <- deltas

  # find best delta
  best_idx <- which.min(mean_errors)
  delta_star <- deltas[best_idx]

  # Refit final NOVELIST covariance using entire base_resid
  novelist_results <- novelist_est(
    resid     = base_resid,
    delta     = delta_star,
    zero_mean = zero_mean
  )

  final_cov_novelist <- novelist_results$cov.novelist
  lambda_star <- novelist_results$lambda

  # Check final positive definiteness
  if (any(eigen(final_cov_novelist)$values <= 1e-8)) {
    stop("The final covariance matrix is not positive definite.")
  }

  # Return results
  return(list(
    delta = delta_star,
    lambda = lambda_star,
    cov.novelist = final_cov_novelist,
    errors_by_delta = mean_errors
  ))
}



#' Reconcile Forecasts using Min Trace approach
#'
#' @param base_forecasts A matrix of base forecasts (dim T by p).
#' @param S A matrix of reconciliation weights (dim p by b, b is number of bottom series).
#' @param W A covariance matrix of h-step base forecast errors (dim p by p).
#'
reconcile_mint <- function(base_forecasts, S, W) {

  R <- t(S)%*%solve(W)
  P <- solve(R%*%S)%*%R

  if (is.vector(base_forecasts)) {
    recon_fc <- S %*% P %*% base_forecasts

  } else if (nrow(S) == ncol(base_forecasts)) {
    # reconcile each row of base forecasts
    recon_fc <- base_forecasts
    for (i in 1:nrow(base_forecasts)) {
      recon_fc[i, ] <- S %*% P %*% base_forecasts[i, ]
    }

  } else {
    stop("This function takes base forecasts as vector or T by p matrix, and S as p by b matrix.")
  }

  return(recon_fc)
}
