# -------------------------------------------------
# ----------- Shrinkage Estimator Tests -----------
# -------------------------------------------------

test_that("shrinkage_est handles normal usage correctly", {
  x <- matrix(rnorm(5 * 20), nrow = 20, ncol = 5)
  out <- shrinkage_est(x)

  # lambda is within [0, 1]
  expect_true(out$lambda >= 0 && out$lambda <= 1)

  # 5x5, positive definite, numeric
  expect_equal(dim(out$cov), c(5,5))
  expect_true(all(is.finite(out$cov)))
  # symmetrical
  expect_equal(out$cov, t(out$cov), tolerance = 1e-14)
  # PD by checking all eigenvalues > 0
  evs <- eigen(out$cov, only.values = TRUE)$values
  expect_true(all(evs > 1e-8))
})


test_that("shrinkage_est with user-specified lambda works", {
  x <- matrix(rnorm(3 * 10), nrow = 10, ncol = 3)
  out <- shrinkage_est(x, lambda = 0.5)

  expect_equal(out$lambda, 0.5)
  expect_equal(dim(out$cov), c(3,3))
})

test_that("shrinkage_est returns stable results for p > n", {
  x <- matrix(rnorm(10*20), nrow = 10, ncol = 20)  # p=20 > n=10

  out <- shrinkage_est(x)
  expect_true(out$lambda >= 0 && out$lambda <= 1)
  # eigen check
  evs <- eigen(out$cov, only.values = TRUE)$values
  expect_true(all(evs > 1e-8))
})

test_that("shrinkage_est stops if lambda out of [0,1]", {
  x <- matrix(rnorm(10), nrow = 5, ncol = 2)
  expect_error(
    shrinkage_est(x, lambda = 1.5)
  )
  expect_error(
    shrinkage_est(x, lambda = -0.2)
  )
})

test_that("check with fabletools function", {
  x <- matrix(rnorm(5 * 20), nrow = 20, ncol = 5)

  shrinkage_est_fable <- function(x) {
    n = nrow(x)
    covm <- crossprod(stats::na.omit(x)) / n # assume x = 0
    tar <- diag(apply(x, 2, purrr::compose(crossprod, stats::na.omit))/n) # SSR/n
    corm <- cov2cor(covm)
    xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
    xs <- xs[stats::complete.cases(xs),]
    v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
    diag(v) <- 0
    corapn <- cov2cor(tar)
    d <- (corm - corapn)^2
    lambda <- sum(v)/sum(d)
    lambda <- max(min(lambda, 1), 0)
    W <- lambda * tar + (1 - lambda) * covm
    return(W)
  }

  out <- shrinkage_est_fable(x)
  myout <- shrinkage_est(x, zero_mean = TRUE)
  expect_equal(out, myout$cov, tolerance = 1e-10)
})

# ------------------------------------------------
# ----------- NOVELIST Estimator Tests -----------
# ------------------------------------------------

test_that("novelist_est basic usage (auto-lambda)", {
  x <- matrix(rnorm(5 * 20), nrow=20, ncol=5)

  out <- novelist_est(x, delta = 0.3)

  # auto-lambda should be in [0, 1]
  expect_true(out$lambda >= 0 && out$lambda <= 1)

  # Cov matrix checks
  expect_equal(dim(out$cov), c(5,5))
  expect_true(all(is.finite(out$cov)))
  expect_equal(out$cov, t(out$cov), tolerance=1e-14)
  evs <- eigen(out$cov, only.values = TRUE)$values
  expect_true(all(evs > 1e-8))
})

test_that("novelist_est works with user-specified lambda", {
  x <- matrix(rnorm(5 * 10), nrow=10, ncol=5)

  out <- novelist_est(x, delta=0.3, lambda=0.6)
  expect_equal(out$lambda, 0.6)
})

test_that("novelist_est warns or stops if delta out of [0,1]", {
  x <- matrix(rnorm(5 * 10), nrow=10, ncol=5)
  expect_error(
    novelist_est(x, delta=-0.5)
  )
  expect_error(
    novelist_est(x, delta=1.1)
  )
})


# TODO: Novelist cannot handle non-positive definite yet
#
# test_that("novelist_est stable with p > n", {
#   x <- matrix(rnorm(10*20), nrow = 10, ncol = 20)  # p=20 > n=10
#   # threshold = 0.5
#   out <- novelist_est(x, delta=0.5)
#   expect_true(out$lambda >= 0 && out$lambda <= 1)
#   evs <- eigen(out$cov, only.values=TRUE)$values
#   expect_true(all(evs > 1e-8))
# })


test_that("check with shrinkage estimator", {
  x <- matrix(rnorm(5 * 20), nrow = 20, ncol = 5)
  out <- shrinkage_est(x)
  out2 <- novelist_est(x, delta = 1)

  expect_equal(out$cov, out2$cov, tolerance = 1e-10)
})
