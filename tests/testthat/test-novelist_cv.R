# test‑novelist‑cv.R
library(testthat)

#####  helper to build a *very* small reconciliation structure  #####
# For unit‑tests we don't need a real hierarchy; an identity S works
# (MinT then returns the base forecasts unchanged, which is easy to assert).
make_identity_S <- function(p) diag(p)

##### 1. Basic functionality and structure ##########################
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
# ##### 3. p > n situation #############################################
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

##### 4. Custom error metric #########################################
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

##### 5. Error handling ##############################################
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
# ##### 6. Positive‑definiteness enforcement inside rolling loop ########
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

