groups <- c(2,3)
k <- length(groups)
p <- sum(groups)
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


# Combine Ablocks into block diagonal A
A <- matrix(0, nrow=p, ncol=p)
idx1 <- 1
for(i in seq_len(k)) {
  gsize <- groups[i]
  idx2 <- idx1 + gsize - 1
  A[idx1:idx2, idx1:idx2] <- Ablocks[[i]]
  idx1 <- idx2 + 1
}





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
