# -------------------------------------------------
# --------------- CV NOVELIST Tests ---------------
# -------------------------------------------------

# For unit‑tests we don't need a real hierarchy; an identity S works
# (MinT then returns the base forecasts unchanged, which is easy to assert).
make_identity_S <- function(p) diag(p)

##### 1. Basic functionality and structure
test_that("novelist_cv returns the right structure and PD covariance", {
  T  <- 30          # length of series
  p  <- 4           # number of series

  y   <- matrix(rnorm(T*p), nrow = T, ncol = p)
  y_hat <- y + matrix(rnorm(T*p, sd = 0.2), nrow = T, ncol = p)

  S_mat <- make_identity_S(p)      # 4x4

  res <- novelist_cv(
    y = y,
    y_hat = y_hat,
    S = S_mat,
    window_size = 10,
    deltas = c(0, 0.25, 0.5, 0.75, 1),
    message = FALSE
  )

  # object structure
  expect_true(is.list(res))
  expect_named(
    res,
    c("delta", "lambda", "cov", "errors_by_delta"),
    ignore.order = TRUE
  )

  # delta chosen must be in candidate set
  expect_true(res$delta %in% c(0, 0.25, 0.5, 0.75, 1))

  # lambda inside [0,1]
  expect_true(res$lambda >= 0 && res$lambda <= 1)

  # covariance matrix properties
  expect_equal(dim(res$cov), c(p, p))
  expect_equal(res$cov, t(res$cov), tolerance = 1e-12)
  expect_true(all(eigen(res$cov)$values > 0))

  # errors_by_delta length matches #deltas, names are correct
  expect_equal(length(res$errors_by_delta), 5)
  expect_equal(as.numeric(names(res$errors_by_delta)), c(0, 0.25, 0.5, 0.75, 1))
})


# TODO: p > window size means non-PD, cannot handle yet
#
# ##### 3. p > n situation 
# test_that("novelist_cv works when p > window_size", {
#   set.seed(3)
#   T <- 25
#   p <- 8
#   y  <- matrix(rnorm(T*p), nrow = T, ncol = p)
#   bf <- y + matrix(rnorm(T*p, sd = 0.3), nrow = T)
#
#   # small window --> p > window (8 > 6)
#   res <- novelist_cv(
#     y, bf, S = make_identity_S(p),
#     window_size = 6,
#     deltas      = c(0, 0.4, 0.8)
#   )
#
#   # basic assertions
#   expect_true(res$lambda >= 0 && res$lambda <= 1)
#   expect_equal(dim(res$cov), c(p, p))
# })

##### 4. Custom error metric 
test_that("novelist_cv works with user‑supplied error_metric (MAE)", {
  y  <- matrix(rnorm(40),  nrow = 10, ncol = 4)
  y_hat <- y + matrix(rnorm(40, sd = 0.15), nrow = 10)

  res <- novelist_cv(
    y, y_hat, S = make_identity_S(4),
    window_size = 4,
    deltas      = c(0, 0.5, 1),
    error_metric = function(a, f) mean(abs(a - f), na.rm = TRUE), # MAE
    message = FALSE
  )

  expect_true(is.numeric(res$errors_by_delta))
  expect_equal(length(res$errors_by_delta), 3)
})

##### 5. Error handling 
test_that("novelist_cv throws informative errors for bad inputs", {
  set.seed(5)
  y_ok  <- matrix(rnorm(30), nrow = 10, ncol = 3)
  y_hat_ok <- y_ok

  S_ok <- make_identity_S(3)

  expect_error(
    novelist_cv(y_ok, y_ok[1:5, ], S_ok, window_size = 5)
  )
  expect_error(
    novelist_cv(y_ok, y_hat_ok, S_ok, window_size = 1)
  )
  expect_error(
    novelist_cv(y_ok, y_hat_ok, S_ok, window_size = 15)
  )
})

# TODO: Novelist cannot handle non-positive definite yet
#
# ##### 6. Positive‑definiteness enforcement inside rolling loop
# test_that("internal PD check inside rolling loop works", {
#   set.seed(6)
#   y  <- matrix(rnorm(50), nrow = 10, ncol = 5)
#   bf <- y
#
#   # Force novelist_est to produce diagonal zero threshold (delta = 0)
#   # but shrinkage ensures PD; we just test it does NOT error
#   expect_silent(
#     novelist_cv(
#       y, bf, S = make_identity_S(5),
#       window_size = 5,
#       deltas = 0        # single candidate
#     )
#   )
# })










# -------------------------------------------------
# ----------- CV NOVELIST with PCA Tests ----------
# -------------------------------------------------

# helper: quick synthetic data
make_demo <- function(T = 40, p = 6, h = 1) {
  set.seed(123)
  # latent 2-factor structure
  F1 <- arima.sim(model = list(ar = 0.6), n = T)
  F2 <- arima.sim(model = list(ar = -0.4), n = T)
  loadings <- matrix(rnorm(p * 2), p, 2)
  eps      <- matrix(rnorm(T * p), T, p) * 0.3
  y  <- cbind(F1, F2) %*% t(loadings) + eps          # “actuals”
  y  <- unclass(y)                                   # plain matrix

  # naïve forecasts = last value + noise
  y_hat <- rbind(NA, y[-T, ]) + matrix(rnorm(T * p, 0, .2), T, p)
  y_hat[1, ] <- y[1, ]                               # no forecast for t=1

  # identity “hierarchy” (for unit test simplicity)
  S <- diag(p)

  list(y = y, y_hat = y_hat, S = S, h = h)
}

# 1. happy path: K = 0,1,2
test_that("novelist_pc_cv returns expected structure", {
  demo <- make_demo(T = 50, p = 8)
  fit <- novelist_pc_cv(
    y  = demo$y,
    y_hat = demo$y_hat,
    Ks = c(0, 1, 2),
    S  = demo$S,
    window_size = 20,
    deltas = c(0, 0.25, 0.5),   # keep test fast
    h = demo$h,
    message = FALSE
  )

  expect_named(fit,
               c("delta", "lambda", "cov", "errors", "ranking_K"))
  expect_length(fit$delta,  3L)          # one per K
  expect_length(fit$lambda, 3L)
  expect_length(fit$cov,    3L)
  # each covariance matrix should be positive-definite
  lapply(fit$cov, function(Sig)
    expect_true(all(eigen(Sig, only.values = TRUE)$values > 0)))
  # errors matrix: rows = |deltas|, cols = |Ks|
  expect_equal(dim(fit$errors), c(3, 3))
})

# 2. dimension mismatch
test_that("mismatched y / y_hat dims throw", {
  demo <- make_demo()
  bad_hat <- demo$y_hat[, -1]           # drop one column
  expect_error(
    novelist_pc_cv(
      y = demo$y, y_hat = bad_hat,
      Ks = 1, S = demo$S,
      window_size = 20, message = FALSE
    ),
    "must match"
  )
})

# 3. invalid window size
test_that("illegal window sizes are rejected", {
  demo <- make_demo()
  expect_error(
    novelist_pc_cv(
      y = demo$y, y_hat = demo$y_hat,
      Ks = 1, S = demo$S,
      window_size = 1,    # too small
      message = FALSE
    ),
    "greater than 1"
  )
  expect_error(
    novelist_pc_cv(
      y = demo$y, y_hat = demo$y_hat,
      Ks = 1, S = demo$S,
      window_size = nrow(demo$y) + 1,   # too large
      message = FALSE
    ),
    "less than nrow"
  )
})

# 4. K = 0, no PCA 
test_that("K = 0 returns no PCA, just NOVELIST CV", {
  demo <- make_demo()
  fit_pc <- novelist_pc_cv(
    y = demo$y, y_hat = demo$y_hat,
    Ks = 0, S = demo$S,
    window_size = 20, deltas = c(0, 0.5),
    h = demo$h, message = FALSE
  )

  fit_pc2 <- novelist_pc_cv(
    y = demo$y, y_hat = demo$y_hat,
    Ks = c(0, 5), S = demo$S,
    window_size = 20, deltas = c(0, 0.5),
    h = demo$h, message = FALSE
  )

  fit <- novelist_cv(
    y = demo$y, y_hat = demo$y_hat,
    S = demo$S,
    window_size = 20, deltas = c(0, 0.5),
    message = FALSE
  )

  expect_true(is.matrix(fit_pc$cov))          # only one covariance matrix 
  expect_equal(length(fit_pc2$cov), 2L)        # list of 2 covariance matrices

  expect_true(all(eigen(fit_pc$cov[[1]], only.values = TRUE)$values > 0))
  expect_true(all(eigen(fit_pc2$cov[[2]], only.values = TRUE)$values > 0))

  expect_equal(fit_pc$cov, fit$cov)
  expect_equal(fit_pc2$cov[[1]], fit$cov)
})

