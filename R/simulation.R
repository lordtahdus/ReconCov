
simulate_bottom_var <- function(
    groups = c(2, 3),      # Vector: size of each group. E.g. c(2, 3) => total 5 series, group1=2, group2=3
    T = 100,               # Number of output time steps
    burnin = 50,           # Burn-in steps (discarded)

    # VAR(1) Coefficients:
    #  - either user-supplied list of length length(groups),
    #    each an (g_i x g_i) matrix,
    #  - or we generate random stable blocks
    Ablocks = NULL,

    # Innovation covariance:
    #  - either user-supplied list of covariance matrices (length(groups)),
    #  - or let the function create them
    Sigblocks = NULL,

    # If we need to generate a block covariance, we specify how:
    corType = c("nonnegative", "mixed"),  # do we allow negative correlation or not
    corRange = c(0.2, 0.7),              # uniform range of correlation magnitudes
    stdevRange = c(sqrt(2), sqrt(6)),    # standard deviations for each series in the group
    # For compound-symmetric approach, each block has all off-diag = correlation,
    # user can choose "mixed" to flip some signs randomly

    random_seed = NULL
) {
  if(!is.null(random_seed)) set.seed(random_seed)

  corType <- match.arg(corType)

  # ----- 1) Basic Dimensions -----
  k <- length(groups)               # number of groups
  p <- sum(groups)                  # total dimension
  # We'll build a block diagonal VAR(1) matrix A of size p x p:

  # ----- 2) Create or Check Ablocks -----
  if(is.null(Ablocks)) {
    # Generate stable block submatrices for each group
    Ablocks <- vector("list", k)
    for(i in seq_along(groups)) {
      gsize <- groups[i]
      # Example: we create an AR coefficient matrix with some random angles/eigenvalues
      # We'll ensure it's stable by bounding spectral radius < 1
      # A simple approach: diagonal of phi's in [0.4, 0.9], off-diag small
      diag_vals <- runif(gsize, min=0.4, max=0.9)
      # off diagonal smaller
      offMat <- matrix(runif(gsize*gsize, min=-0.1, max=0.1), nrow=gsize)
      # Ensure it is stable:
      # We'll do a simple trick: place diag_vals on diagonal, scale the off-diag
      # so that spectral radius is <1
      Ab <- diag(diag_vals) + 0.5 * offMat
      # optional refinement: reduce until stable
      # Not guaranteed stable, so let's do a small loop:
      sr <- max(abs(eigen(Ab)$values))
      while(sr >= 0.99) {
        offMat <- offMat * 0.9
        Ab <- diag(diag_vals) + 0.5 * offMat
        sr <- max(abs(eigen(Ab)$values))
      }
      Ablocks[[i]] <- Ab
    }
  } else {
    # user-supplied, must check dimension
    if(length(Ablocks) != k) {
      stop("Ablocks must be a list of length(groups). Each entry is a g_i x g_i matrix.")
    }
    for(i in seq_along(groups)) {
      if(!all(dim(Ablocks[[i]]) == groups[i])) {
        stop("One of the supplied Ablocks has wrong dimension.")
      }
      # user presumably ensures stability
    }
  }

  # Combine Ablocks into block diagonal A
  A <- matrix(0, nrow=p, ncol=p)
  idx1 <- 1
  for(i in seq_len(k)) {
    gsize <- groups[i]
    idx2 <- idx1 + gsize - 1
    A[idx1:idx2, idx1:idx2] <- Ablocks[[i]]
    idx1 <- idx2 + 1
  }

  # ----- 3) Create or Check Sigblocks (Innovation Covariances) -----
  if(is.null(Sigblocks)) {
    Sigblocks <- vector("list", k)
    for(i in seq_along(groups)) {
      gsize <- groups[i]
      # build a compound-symmetric matrix
      # first pick random correlation in corRange
      r <- runif(1, min=corRange[1], max=corRange[2])
      # might flip sign if corType=="mixed"
      if(corType=="mixed") {
        signflip <- sample(c(-1,1), size=1)
        r <- signflip * r
      }
      # pick stdevs from stdevRange
      sds <- runif(gsize, min=stdevRange[1], max=stdevRange[2])
      # build covariance
      # compound symmetric => Cov(i,i)=sds[i]^2, Cov(i,j)= r*sds[i]*sds[j]
      Sblock <- diag(sds^2)
      for(ii in 1:(gsize-1)) {
        for(jj in (ii+1):gsize) {
          val <- r * sds[ii]*sds[jj]
          Sblock[ii,jj] <- val
          Sblock[jj,ii] <- val
        }
      }
      Sigblocks[[i]] <- Sblock
    }
  } else {
    # user-supplied
    if(length(Sigblocks) != k) {
      stop("Sigblocks must be a list of length(groups). Each entry is g_i x g_i.")
    }
    # no further checks beyond dimension
  }

  # Combine Sigblocks into block diagonal Sig
  Sig <- matrix(0, nrow=p, ncol=p)
  idx1 <- 1
  for(i in seq_len(k)) {
    gsize <- groups[i]
    idx2 <- idx1 + gsize - 1
    Sig[idx1:idx2, idx1:idx2] <- Sigblocks[[i]]
    idx1 <- idx2 + 1
  }

  # ----- 4) Simulate data from the VAR(1) Process -----
  # We do a direct simulation: y_t = A y_{t-1} + eps_t, eps_t ~ N(0, Sig)

  # We'll store entire time series from t=1..(T+burnin)
  Yfull <- matrix(0, nrow=T+burnin, ncol=p)

  # for initialization let's just use zero or small random
  # optionally do a stationarity approach, but let's keep it simple
  # Sample y_1 from stationary dist, or just zero
  # For a robust approach, we do a short burnin anyway.

  # We'll do random initialization:
  Yfull[1,] <- mvrnorm(1, mu=rep(0,p), Sigma=Sig)

  for(t in 2:(T+burnin)) {
    eps_t <- MASS::mvrnorm(1, mu=rep(0,p), Sigma=Sig)
    Yfull[t,] <- A %*% Yfull[t-1,] + eps_t
  }

  # discard the first 'burnin' steps
  Y <- Yfull[(burnin+1):(T+burnin), , drop=FALSE]

  # return results
  list(
    Y = Y,             # final T x p matrix of bottom series
    A = A,             # the block-diagonal VAR(1) coefficient matrix
    Sig = Sig,         # the block-diagonal innovation covariance
    groups = groups,
    Ablocks = Ablocks,
    Sigblocks = Sigblocks,
    burnin = burnin
  )
}
