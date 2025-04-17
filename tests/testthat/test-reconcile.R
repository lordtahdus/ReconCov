test_that("reconcile_mint returns each row unchanged for S=I, W=I with matrix input", {
  T  <- 7
  p  <- 4
  fc_mat <- matrix(rnorm(T*p), nrow = T, ncol = p)
  S  <- diag(p)
  W  <- diag(p)

  out <- reconcile_mint(fc_mat, S, W)

  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(fc_mat))
  expect_equal(out, fc_mat, tolerance = 1e-12)
})

test_that("reconcile_mint enforces coherence in a 2‑level hierarchy", {
  # 4 bottom series, 2 intermediate, 1 top => total m = 7
  S <- rbind(
    # top
    c(1, 1, 1, 1),
    # level 1
    c(1, 1, 0, 0),
    c(0, 0, 1, 1),
    # bottom
    diag(4)
  )
  m <- nrow(S)  # 7
  p <- ncol(S)  # 4
  # As simplest test, use W = I_m
  W <- diag(m)
  base_fc <- seq_len(m)
  out      <- reconcile_mint(base_fc, S, W)
  expect_length(out, m)

  # If W = I, then P = (S'S)^(-1) S', so recon = S P base_fc
  bottoms <- out[4:7]
  lvlA    <- out[2]
  lvlB    <- out[3]
  top     <- out[1]

  expect_equal(lvlA, bottoms[1] + bottoms[2], tolerance = 1e-12)
  expect_equal(lvlB, bottoms[3] + bottoms[4], tolerance = 1e-12)
  expect_equal(top, lvlA + lvlB, tolerance = 1e-12)

  # Now test matrix input: two time points, same base_fc repeated
  fc_mat <- rbind(base_fc, base_fc + 10)
  out_mat <- reconcile_mint(fc_mat, S, W)
  expect_equal(dim(out_mat), c(2, m))
  # rowwise coherence holds
  for(i in 1:2){
    btm <- out_mat[i, 4:7]
    expect_true(abs(out_mat[i,2] - (btm[1] + btm[2])) < 1e-10)
    expect_true(abs(out_mat[i,3] - (btm[3] + btm[4])) < 1e-10)
    expect_true(abs(out_mat[i,1] - (out_mat[i,2] + out_mat[i,3])) < 1e-10)
  }
})


test_that("invalid dimensions", {
  # S is 4x2, but fc is 3x4 → mismatch
  S  <- matrix(1, nrow = 4, ncol = 2)
  W  <- diag(4)
  fc_bad <- matrix(rnorm(3*4), nrow = 3, ncol = 4)
  expect_error(
    reconcile_mint(fc_bad, S, W)
  )
  # fc is a vector too long
  fc_vec_bad <- rnorm(5)
  expect_error(
    reconcile_mint(fc_vec_bad, S, W)
  )
})
