#' Simulate Bottom-Level Series from a Block-Diagonal VAR(1) Model
#'
#' This function simulates a set of bottom-level time series for a hierarchical
#' structure. The series follow a VAR(1) process with block-diagonal coefficients
#' and a specified covariance structure for innovations.
#'
#' @param groups Integer vector specifying the sizes of each bottom-level group, e.g.
#'   \code{c(2,3,4)} means three blocks (2D, 3D, 4D) for a total of \eqn{\sum groups} series.
#' @param T Integer, the number of time steps of final output (excluding burn-in).
#' @param intercept Numeric, the intercept term for the VAR(1) process (default is 0).
#' @param A \eqn{p \times p} matrix of VAR(1) coefficients. If \code{NULL} (default),
#'   a block-diagonal matrix is generated via \code{\link{generate_block_diag}}.
#' @param Sig \eqn{p \times p} covariance matrix of innovations. If \code{NULL}
#'   (default), the function generates a correlation matrix with \code{\link{generate_cor}} and
#'   converts it to covariance. If \code{corType} is "mixed", signs can be flipped.
#' @param burnin Integer, how many initial steps to discard as burn-in (defaults to 50).
#' @param random_seed Optional integer seed for reproducibility. If not \code{NULL}, the
#'   random number generator is set once at the beginning of the function.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{Y}}{A \code{T x p} matrix of simulated bottom-level time series.}
#'     \item{\code{A}}{The final \eqn{p x p} VAR(1) coefficient matrix.}
#'     \item{\code{Sig}}{The final \eqn{p x p} innovation covariance matrix.}
#'     \item{\code{groups}}{The same integer vector specifying group sizes.}
#'   }
#'
#' @export
simulate_bottom_var <- function(
    groups = c(2, 3, 4),
    T = 100,
    intercept = 0,
    A = NULL,
    Sig = NULL,
    burnin = 50,
    random_seed = NULL
) {
  if(!is.null(random_seed)) set.seed(random_seed)

  k <- length(groups)               # number of groups
  p <- sum(groups)                  # total dimension

  # Block matrix of VAR(1) coefficients
  if(is.null(A)) {
    A <- generate_block_diag(
      groups = groups
    )$A
  } else {
    # user-supplied
    if(ncol(A) != p || nrow(A) != p) {
      stop("'A' must be a p x p matrix, where p = sum(groups).")
    }
  }

  # Noise Covariance Matrix
  if (is.null(Sig)) {
    Sigcor <- generate_cor(
      groups = groups
    )
    # Convert correlation to covariance
    Sig <- convert_cor2cov(
      cor = Sigcor
    )
  } else {
    # user-supplied
    if(ncol(Sig) != p || nrow(Sig) != p) {
      stop("'Sig' must be a p x p matrix, matching sum(groups).")
    }
  }

  # Simulate data from the VAR(1) Process
  Yfull <- matrix(0, nrow = T + burnin, ncol = p)
  Yfull[1,] <- MASS::mvrnorm(1, mu=rep(0,p), Sigma=Sig) #random initialise
  eps <- MASS::mvrnorm(T + burnin, mu=rep(0,p), Sigma=Sig) # random noise
  # y_t = A y_{t-1} + eps_t, eps_t ~ N(0, Sig)
  for(t in 2:(T+burnin)) {
    Yfull[t,] <- A %*% Yfull[t-1,] + eps[t,]
  }
  # TODO: optimise

  # discard the first 'burnin' steps
  Y <- Yfull[(burnin+1):(T+burnin), , drop=FALSE]
  Y <- Y + intercept

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
#' @param groups Integer vector specifying the number of series in each group
#' @param rho Numeric vector of length \code{length(groups)} specifying the
#'  base correlation values for each group.
#' @param delta Numeric, the correlation value for off-diagonal blocks (default is 0.15).
#' @param epsilon Numeric, the scaling factor for noise added to the correlation matrix
#'  (default is 0.15).
#' @param eidim Integer, the dimension of the noise space (default is 2).
#' @param random_seed Optional integer for reproducibility.
#' @return A \eqn{p \times p} correlation matrix, where \eqn{p = \sum groups}.
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

  # Constraints on delta and epsilon
  if (delta < 0 || delta >= min(rho)) {
    stop("Delta must be between 0 and minimum of rho parameters.")
  }
  if (epsilon < 0 || epsilon >= (1 - max(rho))) {
    stop("Epsilon must be between 0 and (1 - maximum of rho parameters).")
  }

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
    # TODO: use Matrix::bdiag() instead
  }
  diag(bigcor) <- 1 - epsilon

  eivect <- c()
  for (i in 1:p) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  # TODO: avoid loop
  # E_raw <- matrix(runif(eidim * p, min = -1, max = 1), nrow = eidim)
  # norms <- sqrt(colSums(E_raw^2))
  # eivect <- sweep(E_raw, 2, norms, FUN = "/") * sqrt(epsilon)

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
    stationary = TRUE,
    random_seed = NULL,
    message = TRUE
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

    if (stationary) {
      # Ensure stability by decreasing the scale until spectral radius < 1 or attempts done
      sr <- max(abs(eigen(block)$values))
      attempt <- 0
      # Warning
      if (sr >= 0.95 && message) {
        message("generate_block_diag(): \n",
                "Simulated ", i,"-th block matrix is unstable (not stationary).")
      }
      while(sr >= 0.95) {
        # reduce the off diag portion
        offdiag_mat <- 0.95 * offdiag_mat
        block <- diag(diag_vals) + offdiag_mat
        sr <- max(abs(eigen(block)$values))
        attempt <- attempt + 1
      }
      if (attempt > 0 && message) message("Rescaled by 0.95^", attempt)

      # TODO: ensure max attempts
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
  # check if cor is a correlation matrix
  if (!is.matrix(cor) || nrow(cor) != ncol(cor)) {
    stop("Input 'cor' must be a square matrix.")
  }
  # stays in -1 and 1
  if (any(abs(cor) > 1)) {
    stop("Input 'cor' must be a correlation matrix with values in [-1, 1].")
  }
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
  if (any(eigen(mat_new)$values <= 1e-10)) {
    warning("Flipped matrix is not positive definite.\n")
    # Ensure PD?
    if (ensure_PD) {
      iscorr <- all(abs(diag(mat_new) - 1) < 1e-10)
      mat_new <- Matrix::nearPD(mat_new, corr = iscorr)$mat
      message("Converted to PD matrix using nearPD() algorithm of Higham (2002).\n")
    }
  }

  return(mat_new)
}


#' Combine blocks into a block diagonal matrix
#'
#' TODO: use Matrix::bdiag() |> as.matrix() instead
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


