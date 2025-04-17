# -----------------------------------------------------
# --------------- generate_block_diag() ---------------
# -----------------------------------------------------

test_that("Return list has fields A and blocks with correct dimensions", {
  groups <- c(2, 3, 50)
  p <- sum(groups)
  res <- generate_block_diag(
    groups = groups,
    diag_range = c(0.7, 0.9),
    offdiag_range = c(-0.1, 0.1),
    stationary = TRUE,
    random_seed = 1,
    message = FALSE
  )

  # Total dimension should be 105
  expect_equal(dim(res$A), c(p, p))

  # Blocks list should have 3 matrices with dimensions 2x2, 3x3, and 4x4 respectively
  expect_equal(length(res$blocks), 3)
  expect_equal(dim(res$blocks[[1]]), rep(groups[1], 2))
  expect_equal(dim(res$blocks[[2]]), rep(groups[2], 2))
  expect_equal(dim(res$blocks[[3]]), rep(groups[3], 2))

  # Spectral radius of each block is below 0.99
  for (block in res$blocks) {
    sr <- max(abs(eigen(block)$values))
    expect_lt(sr, 0.99)
  }
})

test_that("A is block-diagonal with each block matching corresponding submatrix", {
  groups <- c(2, 3, 4)
  res <- generate_block_diag(groups = groups, random_seed = 123)
  A <- res$A
  start_idx <- 1
  for (i in seq_along(res$blocks)) {
    block <- res$blocks[[i]]
    gsize <- nrow(block)
    end_idx <- start_idx + gsize - 1
    # The block submatrix should equal the block computed
    submat <- A[start_idx:end_idx, start_idx:end_idx]
    expect_equal(submat, block)

    # Off-block parts in the corresponding rows/cols should be zero:
    if (start_idx > 1) {
      left_block <- A[start_idx:end_idx, 1:(start_idx - 1)]
      expect_true(all(abs(left_block) < 1e-12))
    }
    if (end_idx < ncol(A)) {
      right_block <- A[start_idx:end_idx, (end_idx + 1):ncol(A)]
      expect_true(all(abs(right_block) < 1e-12))
    }
    start_idx <- end_idx + 1
  }
})

test_that("Using offdiag_range = 0 yields diagonal-only blocks", {
  groups <- c(2, 3, 37)
  res <- generate_block_diag(
    groups = groups,
    offdiag_range = 0,
    stationary = FALSE
  )
  for (block in res$blocks) {
    # Off-diagonals must be exactly 0
    expect_equal(block, diag(diag(block)))
  }
})

test_that("Constant diag_range yields blocks with constant diagonal values", {
  groups <- c(69, 2, 3)
  res <- generate_block_diag(
    groups = groups,
    diag_range = c(0.5, 0.5),
    stationary = FALSE
  )
  for (block in res$blocks) {
    expect_equal(diag(block), rep(0.5, nrow(block)))
  }
})

test_that("A single-series block returns a scalar equal to the specified diagonal", {
  res <- generate_block_diag(
    groups = c(1),
    diag_range = c(0.3, 0.3)
  )
  expect_equal(nrow(res$A), 1)
  expect_equal(ncol(res$A), 1)
  expect_equal(res$A[1, 1], 0.3)
})
