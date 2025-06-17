#' Reconcile Forecasts using Min Trace approach
#'
#' @param base_forecasts A matrix of base forecasts (dim T by p).
#' @param S A matrix of reconciliation weights (dim p by b, b is number of bottom series).
#' @param W A covariance matrix of h-step base forecast errors (dim p by p).
#'
#' @export
reconcile_mint <- function(base_forecasts, S, W) {

  R <- t(S)%*%solve(W)
  P <- solve(R%*%S)%*%R

  if (is.vector(base_forecasts)) {
    recon_fc <- base_forecasts
    recon_fc <- S %*% P %*% base_forecasts
    # return as a 1 by p vector to ensure consistency with matrix form
    recon_fc <- t(recon_fc)

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
