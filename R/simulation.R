
simulate_bottom_var <- function(
    groups = c(2, 3, 4),   # Vector: size of each group. E.g. c(2, 3) => total 5 series, group1=2, group2=3
    T = 100,               # Number of output time steps
    A = NULL,
    Sig = NULL,

    diag_range_A = c(0.4, 0.9),
    offdiag_range_A = c(-0.1, 0.1),

    corType = c("nonnegative", "mixed"),  # do we allow negative correlation or not
    cor_range = c(0.2, 0.7),              # uniform range of correlation magnitudes
    stdev_range = c(sqrt(2), sqrt(6)),    # standard deviations for each series in the group
    burnin = 50,           # Burn-in steps (discarded)
    random_seed = NULL
) {
  if(!is.null(random_seed)) set.seed(random_seed)

  corType <- match.arg(corType)

  k <- length(groups)               # number of groups
  p <- sum(groups)                  # total dimension

  # ----- Block matrix of VAR(1) coefficients -----
  # If no A, call generate_block_diag() function
  if(is.null(A)) {
    A <- generate_block_diag(
      groups = groups,
      diag_range = diag_range_A,
      offdiag_range = offdiag_range_A,
      random_seed = random_seed
    )$A
  } else {
    # user-supplied
    if(ncol(A) != p || nrow(A) != p) {
      stop("Matrix A must have dimension of p x p.")
    }
  }

  # ----- Noise Covariance Matrix -----
  # If no Sig, call generate_cor() function
  if (is.null(Sig)) {
    Sigcor <- generate_cor(
      groups = groups,
      rho = runif(k, min = cor_range[1], max = cor_range[2]),
      delta = min(cor_range) * 0.9,
      epsilon = 0.99 - max(cor_range),
      eidim = 2
    )
    # Convert correlation to covariance
    Sig <- convert_cor2cov(
      cor = Sigcor,
      stdevs = runif(p, min = stdev_range[1], max = stdev_range[2])
    )

    if (corType == "mixed") {
      # convert to mixed sign matrix
      Sig <- flip_signs_mat(Sig, ensure_PD = TRUE)
    }

  } else {
    # user-supplied
    if(ncol(Sig) != p || nrow(Sig) != p) {
      stop("Error cov martix must have dimension of p x p.")
    }
  }

  # ----- Simulate data from the VAR(1) Process -----
  # y_t = A y_{t-1} + eps_t, eps_t ~ N(0, Sig)

  Yfull <- matrix(0, nrow = T + burnin, ncol = p)
  # Random initialise
  Yfull[1,] <- MASS::mvrnorm(1, mu=rep(0,p), Sigma=Sig)

  for(t in 2:(T+burnin)) {
    eps_t <- MASS::mvrnorm(1, mu=rep(0,p), Sigma=Sig)
    Yfull[t,] <- A %*% Yfull[t-1,] + eps_t
  }

  # discard the first 'burnin' steps
  Y <- Yfull[(burnin+1):(T+burnin), , drop=FALSE]

  return(list(
    Y = Y,
    A = A,
    Sig = Sig,
    groups = groups
  ))
}


#' Generate Correlation Structure using Hardin et al. (2013) algorithm 1
#'
#' Start with constant correlations within each group (rhos) and between
#' each group (delta). Then add random noise to all elements of the matrix, by
#' scaling epsilon by outer cross products of random unit vectors.
#'
#' @references
#' Hardin, J., Garcia, S. R., & Golan, D. (2013). A method for generating realistic correlation matrices. The Annals of Applied Statistics, 7(3), 1733–1762. https://www.jstor.org/stable/23566492
#'
#' @export
generate_cor <- function (
    groups = c(2, 3, 4),
    rho = runif(length(groups), 0.2, 0.7),
    delta = min(rho) * 0.8,
    epsilon = (1 - max(rho)) * 0.8,
    eidim = 2,
    random_seed = NULL
) {
  if(!is.null(random_seed)) set.seed(random_seed)

  k <- length(groups)               # number of groups
  p <- sum(groups)                  # total dimension

  bigcor <- matrix(rep(delta, p * p), ncol = p)
  for (i in 1:k) {
    cor <- matrix(rep(rho[i], groups[i] * groups[i]), ncol = groups[i])
    if (i == 1) {bigcor[1:groups[1], 1:groups[1]] <- cor}
    if (i != 1) {
      bigcor[
        (sum(groups[1:(i - 1)]) + 1):sum(groups[1:i]),
        (sum(groups[1:(i - 1)]) + 1):sum(groups[1:i])
      ] <- cor
    }
  }
  diag(bigcor) <- 1 - epsilon

  eivect <- c()
  for (i in 1:p) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  bigE <- t(eivect) %*% eivect
  finalcor <- bigcor + bigE
  return(finalcor)
}


#' Generate Block-Diagonal VAR(1) Coefficient Matrix
#'
#' Constructs a block-diagonal matrix for a VAR(1) process, with each block
#' corresponding to a group of series. Blocks are either purely diagonal or
#' contain random off-diagonal entries scaled to ensure stationarity.
#'
#' @param groups Integer vector specifying number of series in each group
#'   (e.g., \code{groups = c(2, 3, 4)} for three blocks of sizes 2, 3, and 4).
#' @param diag_range Numeric vector length 2 specifying range of diagonal entries.
#'   Default is \code{c(0.2, 0.9)}.
#' @param offdiag_range Numeric vector of length 2 specifying the range for
#'   off-diagonal entries, or 0 for diagonal only. Default is \code{c(-0.1, 0.1)}.
#' @param random_seed Optional integer for reproducibility.
#' @return A list with:
#' \describe{
#'   \item{\code{A}}{Full block-diagonal VAR(1) matrix.}
#'   \item{\code{blocks}}{List of individual block matrices.}
#' }
#'
#' @export
generate_block_diag <- function(
    groups = c(2, 3, 4),
    diag_range = c(0.2, 0.9),
    offdiag_range = c(-0.1, 0.1),
    random_seed = NULL
) {
  # Allows user to fix the random seed for reproducibility
  if(!is.null(random_seed)) set.seed(random_seed)

  k <- length(groups)
  p <- sum(groups)   # total dimension

  blocks <- vector("list", k)

  # Loop over each group i
  for(i in seq_len(k)) {
    gsize <- groups[i]
    diag_vals <- runif(gsize, min = diag_range[1], max = diag_range[2])

    # If offdiag_range is 0 or NULL, skip off-diagonal generation
    if (is.null(offdiag_range) || all(offdiag_range == 0)) {
      offdiag_mat <- matrix(0, nrow = gsize, ncol = gsize)
    } else {
      offdiag_mat <- matrix(
        runif(gsize*gsize, min=offdiag_range[1], max=offdiag_range[2]),
        nrow = gsize
      )
      diag(offdiag_mat) <- 0  # zero out diagonal
    }
    # Combine diag & offdiag
    block <- diag(diag_vals, nrow = gsize) + offdiag_mat

    # Ensure stability by decreasing the scale until spectral radius < 1 or attempts done
    sr <- max(abs(eigen(block)$values))
    # Warning
    if (sr >= 0.99) {
      cat("Simulated block matrix is unstable (not stationary).",
          "Attempt to rescale of off-diagonal elements.")
    }
    while(sr >= 0.99) {
      # reduce the off diag portion
      offdiag_mat <- 0.95 * offdiag_mat
      block <- diag(diag_vals) + offdiag_mat
      sr <- max(abs(eigen(block)$values))
    }
    blocks[[i]] <- block
  }
  # Combine blocks into single block diagonal A
  A <- combine_blocks(blocks)

  return(list(
    A = A,             # full p x p matrix
    blocks = blocks    # list of each block
  ))
}


#' Convert correlation to covariance matrix
#'
convert_cor2cov <- function(
    cor,
    stdevs = runif(nrow(cor), sqrt(2), sqrt(6))
) {
  p <- nrow(cor)
  cov <- diag(stdevs) %*% cor %*% diag(stdevs)
  return(cov)
}


#' Flip totally positive matrix to mixed-sign matrix
#'
#' Randomly flips signs of the elements with probability \code{flip_prob}.
#'
flip_signs_mat <- function(mat, flip_prob = 0.5, ensure_PD = TRUE) {
  # Check that mat is a square matrix
  if(nrow(mat) != ncol(mat)) {
    stop("Input matrix must be square.")
  }
  n <- nrow(mat)
  mat_new <- mat
  # Indices for upper & lower triangle (excluding diagonal)
  upper_idx <- which(upper.tri(mat_new))
  lower_idx <- which(lower.tri(mat_new))

  # Generate random multipliers (-1 with probability flip_prob, otherwise +1)
  flip_factors <- ifelse(runif(length(upper_idx)) < flip_prob, -1, 1)

  mat_new[upper_idx] <- mat_new[upper_idx] * flip_factors
  mat_new[lower_idx] <- t(mat_new)[lower_idx]

  # Check PD
  if (any(eigen(mat_new)$values <= 1e-8)) {
    cat("Flipped matrix is not positive definite.\n")
    # Ensure PD?
    if (ensure_PD) {
      iscorr <- all(abs(diag(mat_new) - 1) < 1e-8)
      mat_new <- Matrix::nearPD(mat_new, corr = iscorr)$mat
      cat("Converted to PD matrix using nearPD() algorithm of Higham (2002).\n")
    }
  }

  return(mat_new)
}


#' Combine blocks into a block diagonal matrix
#'
combine_blocks <- function(blocks) {
  k <- length(blocks)
  p <- sum(sapply(blocks, nrow))   # total dimension

  A <- matrix(0, nrow=p, ncol=p)
  idx1 <- 1
  for(i in seq_len(k)) {
    gsize <- nrow(blocks[[i]])
    idx2 <- idx1 + gsize - 1
    A[idx1:idx2, idx1:idx2] <- blocks[[i]]
    idx1 <- idx2 + 1
  }
  return(A)
}


