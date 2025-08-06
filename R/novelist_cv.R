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
#' @param h The forecast horizon from 1 to h (default is 1).
#' @param reconcile_forecast .....
#' @param zero_mean Logical, whether to treat the residuals as zero mean in the covariance.
#' @param error_metric An error measure function for given actual and reconciled forecasts.
#' @param ensure_PD Logical, whether to ensure the covariance matrix is positive definite.
#'
#' @return A list containing:
#'   \item{delta}{optimal threshold}
#'   \item{lambda}{optimal shrinkage intensity}
#'   \item{cov}{optimal NOVELIST covariance matrix}
#'   \item{errors_by_delta}{vector of average validation error for each delta}
#'
#' @details
#' This approach loops through each possible validation segment in a rolling manner:
#'   for i in window_size:(T-1), use resid[(i-window_size+1):i, ] as training data
#'   then produce a reconciled forecast for time i+1 and measure the error.
#' The final best delta is the one that yields the lowest error measure.
#' 
#' @import Matrix
#' @export
novelist_cv <- function(
    y,
    y_hat,
    S,
    window_size,
    deltas = seq(0, 1, by = 0.05),
    h = 1,
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

  # Check possible semi-positive-definite covariance matrix
  if (window_size < p + 1 & !ensure_PD) {
    stop("Sample size condition per window not satisfied: window_size must be at least p + 1 to ensure a positive-definite covariance matrix. Please use ensure_PD = TRUE in order to reconcile.")
  }

  # Calculate base residuals
  resid <- y - y_hat
  # Store errors
  cv_errors <- matrix(NA, nrow = T - window_size, ncol = length(deltas))

  # Rolling over each possible validation step
  # i means the "last index" of the training set is i.
  # Then the "validation" point is i+1 to i+h.
  for (i in window_size:(T - h)) {

    # Training residuals from (i-window_size+1) to i
    train_resid <- resid[(i - window_size + 1):i, , drop=FALSE] # drop=F keep matrix structure

    # The actual data at time i+1 to i+h
    actual_next <- y[(i+1) : (i+h), ]
    fitted_next <- y_hat[(i+1) : (i+h), ]

    # Now loop over candidate threshold delta
    for(delta in deltas){

      # Estimate a NOVELIST covariance using the train_resid
      cov_novelist <- novelist_est(
        resid     = train_resid,
        delta     = delta,
        zero_mean = zero_mean,
        ensure_PD = ensure_PD
      )$cov

      # if (any(eigen(cov_novelist)$values <= 1e-12)) {
      #   # stop(sort(eigen(cov_novelist)$values)[1])
      #   stop("The covariance matrix is not positive definite, cannot reconcile. Try ensure_PD = T")
      # }

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
    # if (message) {
    #   if (idx_i %% 5 == 0 || idx_i == (T - window_size)) {
    #     cat("Iteration: ", idx_i, " / ", (T - window_size), "\n")
    #   }
    # }
  }
  # cat("\n")

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



# novelist_cv_test <- function(
#     y,
#     y_hat,
#     S,
#     window_size,
#     deltas = seq(0, 1, by = 0.05),
#     h = 1,
#     # reconcile_forecast = NULL,
#     zero_mean = TRUE,
#     error_metric = function(actual, fc) mean((actual - fc)^2, na.rm = TRUE),
#     ensure_PD = TRUE,
#     message = TRUE
# ) {
#   T <- nrow(y)
#   p <- ncol(y)

#   # Preliminary checks
#   if(any(dim(y) != dim(y_hat))) {
#     stop("Dimensions of 'y' and 'y_hat' must match.")
#   }
#   if(T <= window_size || window_size <= 1 || (window_size %% 1 != 0)) {
#     stop("Window size must be integer, greater than 1 and less than nrow.")
#   }

#   # Calculate base residuals
#   resid <- y - y_hat
#   # Store errors
#   cv_errors <- matrix(NA, nrow = T - window_size, ncol = length(deltas))

#   # Rolling over each possible validation step
#   # i means the "last index" of the training set is i.
#   # Then the "validation" point is i+1 to i+h.
#   for (i in window_size:(T - h)) {

#     # Training residuals from (i-window_size+1) to i
#     train_resid <- resid[(i - window_size + 1):i, , drop=FALSE] # drop=F keep matrix structure

#     # The actual data at time i+1 to i+h
#     actual_next <- y[(i+1) : (i+h), ]
#     fitted_next <- y_hat[(i+1) : (i+h), ]

#     # Now loop over candidate threshold delta
#     for(d in seq_along(deltas)){
#       delta <- deltas[d]
#       # Estimate a NOVELIST covariance using the train_resid
#       cov_novelist <- novelist_est(
#         resid     = train_resid,
#         delta     = delta,
#         zero_mean = zero_mean,
#         ensure_PD = ensure_PD
#       )$cov

#       if (any(eigen(cov_novelist)$values <= 1e-12)) {
#         # stop(sort(eigen(cov_novelist)$values)[1])
#         stop("The covariance matrix is not positive definite, cannot reconcile. Try ensure_PD = T")
#       }

#       recon_fc <- fitted_next
#       recon_fc <- reconcile_mint(
#         base_forecasts = fitted_next,
#         S = S,
#         W = cov_novelist
#       )

#       # Compute error measure for time i+1:
#       err_val <- error_metric(actual_next, recon_fc)

#       # store
#       idx_i <- i - window_size + 1
#       cv_errors[idx_i, d] <- err_val
#     } # end loop over deltas

#     # TODO: Refine Progress Printing
#     if (message) {
#       if (idx_i %% 5 == 0 || idx_i == (T - window_size)) {
#         cat("Iteration: ", idx_i, " / ", (T - window_size), "\n")
#       }
#     }
#   }
#   cat("\n")

#   # Now compute average error across all rolling iterations for each delta
#   mean_errors <- colMeans(cv_errors, na.rm = TRUE)
#   names(mean_errors) <- deltas

#   # find best delta
#   best_idx <- which.min(mean_errors)
#   delta_star <- deltas[best_idx]

#   # Refit final NOVELIST covariance using entire resid
#   novelist_results <- novelist_est(
#     resid     = resid,
#     delta     = delta_star,
#     zero_mean = zero_mean
#   )

#   final_cov_novelist <- novelist_results$cov
#   lambda_star <- novelist_results$lambda

#   # Check final positive definiteness
#   if (any(eigen(final_cov_novelist)$values <= 1e-8)) {
#     stop("The final covariance matrix is not positive definite.")
#   }

#   # Return results
#   return(list(
#     delta = delta_star,
#     lambda = lambda_star,
#     cov = final_cov_novelist,
#     errors_by_delta = mean_errors
#   ))
# }


#' Rolling Cross-Validation for NOVELIST with PCA Threshold Selection
#' 
#' This function performs a rolling window approach to select the threshold 'delta'.
#' For each window, it applies PCA to the residuals to extract K largest factors,
#' then applies NOVELIST covariance estimation on the remaindar part.
#' 
#' It can take in multiple candidate K values, and, for each K, returns the best 
#' threshold 'delta' based on minimising some measure of forecast error on validation windows.
#' 
#' @param y A matrix of actual values (dim T by p),
#' @param y_hat A matrix of fitted values (dim T by p).
#' @param Ks A numeric vector of candidate number of factors to extract from PCA.
#'           If there is 0, it will not apply PCA and use the full residuals.
#'           Default is 1.
#' @param S A matrix of reconciliation structure (dim p by b; b is number of bottom series).
#' @param window_size The length of the rolling window n.
#' @param deltas A numeric vector of candidate threshold values in [0,1].
#' @param h The forecast horizon from 1 to h (default is 1).
#' @param reconcile_forecast .....
#' @param zero_mean Logical, whether to treat the residuals as zero mean in the covariance.
#' @param error_metric An error measure function for given actual and reconciled forecasts.
#' @param ensure_PD Logical, whether to ensure the covariance matrix is positive definite.
#'
#' @return A list containing:
#'  \item{delta}{vector of optimal threshold for each K}
#'  \item{lambda}{vector of optimal shrinkage intensity for each K}
#'  \item{cov}{list of optimal NOVELIST covariance matrices for each K}
#'  \item{errors}{matrix of average validation error for each delta and K}
#'  \item{ranking_K}{vector of minimum mean errors for each K}
#' 
#' @details
#' Write soon
#' 
#' @import Matrix
#' @export
novelist_pc_cv <- function(
    y,
    y_hat,
    Ks = 1,
    S,
    window_size,
    deltas = seq(0, 1, by = 0.05),
    h = 1,
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

  # Check possible semi-positive-definite covariance matrix
  if (window_size < p + 1 & !ensure_PD) {
    stop("Sample size condition per window not satisfied: window_size must be at least p + 1 to ensure a positive-definite covariance matrix. Please use ensure_PD = TRUE in order to reconcile.")
  }

  # Calculate base residuals
  resid <- y - y_hat
  # Store errors, a list of matrices for each K
  cv_errors <- lapply(
    Ks,
    function(i) matrix(NA, nrow = T - window_size, ncol = length(deltas))
  )
  names(cv_errors) <- paste0("K", Ks)

  # Rolling over each possible validation step
  # i means the "last index" of the training set is i.
  # Then the "validation" point is i+1 to i+h.
  for (i in window_size:(T - h)) {

    # Training residuals from (i-window_size+1) to i
    train_resid <- resid[(i - window_size + 1):i, , drop=FALSE] # drop=F keep matrix structure

    # The actual data at time i+1 to i+h
    actual_next <- y[(i+1) : (i+h), ]
    fitted_next <- y_hat[(i+1) : (i+h), ]

    # Now loop over candidate threshold delta
    for(delta in deltas){

      # Estimate a NOVELIST covariance using the train_resid
      # for each K
      cov_novelist_list <- lapply(
        Ks,
        function(K) {
          # If K is 0, use the full residuals without PCA
          if (K == 0) {
            novelist_est(
              resid     = train_resid,
              delta     = delta,
              zero_mean = zero_mean,
              ensure_PD = ensure_PD
            )$cov
          } else {
            # Otherwise, apply PCA and then NOVELIST
            novelist_pc_est(
              resid     = train_resid,
              K         = K,
              delta     = delta,
              zero_mean = zero_mean,
              ensure_PD = ensure_PD
            )$cov
          }
        }
      )
      
      # TODO: Check PD
      # if (any(eigen(cov_novelist)$values <= 1e-12)) {
      #   stop("The covariance matrix is not positive definite, cannot reconcile. Try ensure_PD = T")
      # }

      recon_fc_list <- lapply(
        cov_novelist_list,
        function(W) {
          reconcile_mint(
            base_forecasts = fitted_next,
            S = S,
            W = W
          )
        }
      )

      # Compute error measure for time i+1 for each K
      err_vals <- sapply(
        recon_fc_list,
        function(recon_fc) error_metric(actual_next, recon_fc)
      )

      # store
      idx_i <- i - window_size + 1
      idx_d <- which(deltas == delta)
      for (k in seq_along(Ks)) {
        idx_k <- which(Ks == Ks[k])
        cv_errors[[idx_k]][idx_i, idx_d] <- err_vals[k]
      }

    } # end loop over deltas

    # TODO: Refine Progress Printing
    # if (message) {
    #   if (idx_i %% 5 == 0 || idx_i == (T - window_size)) {
    #     cat("Iteration: ", idx_i, " / ", (T - window_size), "\n")
    #   }
    # }
  }
  # cat("\n")

  # Now compute average error across all rolling iterations
  # for each K (col) and delta (row)
  mean_errors <- sapply(
    cv_errors,
    function(err_mat) colMeans(err_mat, na.rm = TRUE)
  ) # return matrix, or vector if only 1 K and/or 1 delta
  
  if (is.vector(mean_errors)) {
    mean_errors <- matrix(mean_errors, nrow = 1)
  }
  row.names(mean_errors) <- deltas

  # Find best delta for each K (col)
  delta_stars <- apply(
    mean_errors,
    2,
    function(errs) {
      best_idx <- which.min(errs)
      deltas[best_idx]
    }
  )

  # Refit final NOVELIST covariance for each K, using entire resid
  novelist_results <- lapply(
    Ks,
    function(K) {
      if (K == 0) {
        return(novelist_est(
          resid     = resid,
          delta     = delta_stars[which(Ks == K)],
          zero_mean = zero_mean,
          ensure_PD = ensure_PD
        ))
      }
      novelist_pc_est(
        resid     = resid,
        K         = K,
        delta     = delta_stars[which(Ks == K)],
        zero_mean = zero_mean,
        ensure_PD = ensure_PD
      )
    }
  )
  names(novelist_results) <- paste0("K", Ks)

  final_cov_novelist <- lapply(
    novelist_results,
    function(res) res$cov
  )
  lambda_star <- sapply(
    novelist_results,
    function(res) res$lambda
  )

  # Ensure all covariance matrices are positive definite
  for (cov in final_cov_novelist) {
    if (any(eigen(cov)$values <= 1e-8)) {
      stop("One of the final covariance matrices is not positive definite.")
    }
  }

  # find lowest mean errors for each K (column)
  ranking_K <- sort(
    apply(mean_errors, 2, min)
  )

  # If only one K, convert final_cov_novelist to a vector
  if (length(Ks) == 1) {
    final_cov_novelist <- final_cov_novelist[[1]]
  }
  # Return results
  return(list(
    delta = delta_stars,
    lambda = lambda_star,
    cov = final_cov_novelist,
    errors = mean_errors,
    ranking_K = ranking_K
  ))
}