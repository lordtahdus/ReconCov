library(MASS)
# library(matrixcalc)
library(Matrix)
library(tidyr)
library(ggplot2)

library(fabletools)
library(fable)
library(feasts)
library(tsibble)
library(dplyr)

library(devtools)
load_all()


# Parameters -----------------------------------

# groups <- c(2,2)
# groups <- c(4,4,4,4)
# groups <- c(6,6,6,6,6,6)
groups <- c(50,50)

T <- 304
h <- 4
Tsplit <- T - h

structure <- list(
  groups,
  as.list(seq(1,length(groups))),
  list(c(1,2))
)
# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   list(c(1,2,3,4))
# )
# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   # list(c(1,2,3), c(4,5,6)),
#   # list(c(1,2))
#   list(1:6)
# )
structure <- list(
  rep(10,10),
  as.list(1:10),
  list(1:4, 5:6, 7:10),
  list(1:3)
)

(S <- construct_S(
  structure = structure,
  sparse = FALSE,
  ascending = FALSE
))
order_S <- rownames(S)

# ranges for coefs in VAR
diag_range <- c(-0.45, 0.45)
offdiag_range <- c(-0.2, 0.2)

## VAR(1) block -------------------------
A <- generate_block_diag(
  groups = groups,
  diag_range = diag_range,
  offdiag_range = offdiag_range,
  stationary = TRUE,
)$A

plot_heatmap(A, TRUE)
A <- edit(A) ;colnames(A) <- NULL # edit manually

# check stationary
for (block in seq_along(groups)) {
  size <- groups[block]
  index <- seq(sum(groups[1:(block-1)]) + 1, sum(groups[1:block]))
  print(
    any(abs(A[block, block]) > 0.9)
  )
}


## Sigma ---------------------------------
rho <- runif(length(groups), 0.5, 0.8)
Sigma <- generate_cor(
  groups = groups,
  rho = rho,
  delta = min(rho) * 0.5,
  # delta = 0.15,
  epsilon = (1-max(rho))*0.5,
  # epsilon = 0.15,
  eidim = length(groups)
)

plot_heatmap(Sigma, TRUE)

Sigma <- edit(Sigma) ; colnames(Sigma) <- NULL # edit manually
Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
any((eigen(Sigma)$values) < 1e-8)

# convert to cov using random sd
Sigma <- convert_cor2cov(
  Sigma,
  stdevs = runif(nrow(Sigma), 1*sqrt(2), 2*sqrt(3))
)
# flip signs
V <- diag(x = sample(
  c(1,1),
  size = sum(groups), replace = TRUE,
  prob= c(0.5, 0.5)
))
Sigma <- V %*% Sigma %*% V

plot_heatmap(Sigma %>% cov2cor(), TRUE)

### modify a range of values---------
R <- cov2cor(Sigma)
# scale down by x
lower <- 0.66
upper <- 0.69
R[lower < abs(R) & abs(R) < upper] <-
  R[lower < abs(R) & abs(R) < upper] * 0.05   # scale x

plot_heatmap(R)
any((eigen(R)$values) < 1e-8)
R <- nearPD(R)$mat %>% as.matrix()
D <- diag(sqrt(diag(Sigma)))
Sigma <- D %*% R %*% D

### formulate as NOVELIST ------------
delta <- 0.5
lambda <- 0.2
R <- cov2cor(Sigma)
R_thresh <- sign(R) * pmax(abs(R) - delta, 0)
diag(R_thresh) <- 1
R_novelist <- lambda * R_thresh + (1 - lambda) * R
D <- diag(sqrt(diag(Sigma)))
Sigma <- D %*% R_novelist %*% D


# temporary save
params <- list(
  groups = groups,
  S      = S,
  A      = A,
  Sigma  = Sigma
)
saveRDS(params, "sim/temp_params.rds")

# Function -----------------------------------
run <- function(A = NULL, Sigma = NULL, message = F) {

  # # # # # #
  # generate bottom-up series and transforming data
  bottom <- simulate_bottom_var(groups, T, intercept = 100, A=A, Sig=Sigma)
  hts_mat <- bottom$Y %*% t(S)

  hts <- hts_mat %>%
    as_tibble() %>%
    mutate(time = seq(1, nrow(hts_mat))) %>%
    # select(time, everything()) %>%
    as_tsibble(index = time) %>%
    pivot_longer(
      cols = -time,
      names_to = "series",
      values_to = "value"
    )

  # # # # # #
  # Fit and base forecasts
  fit <- hts %>%
    filter(time <= Tsplit) %>%
    model(
      arima = ARIMA(value)
    )
  fc <- fit |>
    forecast(h = h) |>
    as_fable(response = "value", distribution = value)
  fc <- fc %>%
    mutate(h = time - Tsplit)

  # Get data
  y <- hts_mat[1:Tsplit, order_S]
  actual <- hts_mat[(Tsplit + 1):T, order_S]

  y_hat <- fit %>%
    augment() %>%
    select(series, .fitted) %>%
    pivot_wider(names_from = series, values_from = .fitted, names_sort = FALSE) %>%
    as_tibble() %>%
    select(-time) %>%
    as.matrix()

  base_fc <- fc %>%
    as_tibble() %>%
    select(series, .mean, time) %>%
    pivot_wider(names_from = series, values_from = .mean) %>%
    select(-time) %>%
    as.matrix()

  y_hat <- y_hat[, order_S]
  base_fc <- base_fc[, order_S]

  # # # # # #
  # Get covariance estimates
  W_shr <- shrinkage_est(
    y - y_hat
  )
  window <- round(Tsplit * 0.7)
  # W_n <- novelist_cv(
  #   y,
  #   y_hat,
  #   S,
  #   window = window,
  #   deltas = seq(0, 1, by = 0.05),
  #   ensure_PD = TRUE,
  #   message = message
  # )
  # C++ version
  W_n <- novelist_cv_cpp(
    y,
    y_hat,
    S,
    window = window,
    deltas = seq(0, 1, by = 0.05),
    ensure_PD = TRUE
  )

  # # # # # #
  # Reconcile
  recon_mint_shr <- reconcile_mint(base_fc, S, W_shr$cov)
  recon_mint_n <- reconcile_mint(base_fc, S, W_n$cov)

  sample_cov <- compute_cov_matrix(y - y_hat, zero_mean = T)
  if (any(eigen(sample_cov)$values < 1e-10)) {
    cat("Sample covariance for mint_sample is singular, using nearPD\n")
    sample_cov <- nearPD(sample_cov)$mat
  }
  recon_mint_sample <- reconcile_mint(base_fc, S, sample_cov)

  # Sigma_true <- S %*% Sigma %*% t(S)
  # Sigma_true <- nearPD(Sigma_true)$mat # ensure positive-definite
  # recon_mint_true <- reconcile_mint(base_fc, S, Sigma_true)

  recon_ols <- reconcile_mint(base_fc, S, diag(rep(1, nrow(S)))) # identity matrix

  # # # # # #
  # Return
  SSE <- list(
    base = ((actual - base_fc)^2),
    ols = ((actual - recon_ols)^2),
    mint_shr = ((actual - recon_mint_shr)^2),
    mint_n = ((actual - recon_mint_n)^2),
    mint_sample = ((actual - recon_mint_sample)^2)
    # mint_true = ((actual - recon_mint_true)^2)
  )

  list(
    SSE = SSE,
    W_shr = W_shr$lambda,
    W_n = c(W_n$lambda, W_n$delta),
    W1_hat = compute_cov_matrix(y - y_hat, zero_mean = TRUE)
  )
}


# ---- PARALELL ----------------------------------
library(future.apply)
library(progressr)

handlers(global = TRUE) # Setup progress bar handler
handlers("txtprogressbar")  # or "progress" for a fancier bar

plan(multisession, workers = parallel::detectCores() - 1)

M <- 50

# PARALLEL
# res_list <- future_lapply(seq_len(M), function(i) run(), future.seed=TRUE)

# parallel with progress bar
with_progress({
  p <- progressor(along = 1:M)  # auto sets steps = length

  set.seed(2)
  res_list <- future_lapply(
    X = 1:M,
    FUN = function(i) {
      # this guarantees the bar advances even on error
      on.exit(p(sprintf("Sim %d", i)), add = TRUE)

      ## run the simulation but swallow any error
      tryCatch(
        run(A, Sigma, message = FALSE),

        ## put the error into the return value instead of stopping
        error = function(e) {
          structure(list(message = e$message,
                         call    = e$call,
                         sim_id  = i),
                    class = "sim_error")
        }
      )
    },
    future.seed = TRUE
  )
})

# remove any error simulation
cat("Any error in sim:", any(sapply(res_list, inherits, "sim_error")))
res_list <- res_list[!sapply(res_list, inherits, "sim_error")]

model_names <- c("base", "ols", "mint_shr", "mint_n", "mint_sample")
SSE_cum <- setNames(
  lapply(model_names, function(name) {
    matrix(0, h, length(order_S), dimnames = list(1:h, order_S))
  }),
  model_names
)

# check orders of results
res_list[[1]]$SSE |> names() == SSE_cum |> names()

# wrangle
SSE_cum <- Reduce(function(acc, res) Map(`+`, acc, res$SSE),
                      res_list, init = SSE_cum)
W_shr_store <- sapply(res_list, `[[`, "W_shr")
W_n_store <- t(sapply(res_list, `[[`, "W_n")) ; colnames(W_n_store) <- c("lambda", "delta")

MSE <- lapply(SSE_cum, function(mat) mat / M)

plan(sequential) # Reset to sequential


# Warning message:
# In sqrt(diag(best$var.coef)) : NaNs produced
#
# This happens in fitting ARIMA, when the sum of AR coefs is very
# close to 1, causing difficulties in computing the standard errors.


## benchmark ---------------------
# W_shr_store <- numeric(M)
# W_n_store   <- matrix(0, M, 2,
#                       dimnames = list(NULL, c("lambda", "delta")))
# run_sequential <- function() {
#   for(i in seq_len(M)) {
#     cat("Iteration ", i, "\n")
#
#     res <- run(A, Sigma)
#     # add up the SSE matrices
#     SSE_cum <- Map(`+`, SSE_cum, res$SSE)
#     # store the parameters
#     W_shr_store[i]         <- res$W_shr
#     W_n_store[i, ]         <- res$W_n
#   }
# }
#
# bench::mark(
#   future_lapply(seq_len(M), function(i) run(A, Sigma), future.seed=TRUE),
#   run_sequential(),
#   iterations = 1,
#   check = FALSE
# )



# Save --------------------------------


# Combine all your outputs into a named list
results <- list(
  groups = groups,
  S      = S,
  A      = A,
  Sigma  = Sigma,
  MSE    = MSE,    # or averaged MSE
  W_shr  = W_shr_store,       # can be a list of matrices or one matrix
  W_n    = W_n_store          # same
)
error_list <- purrr::map(res_list, "SSE")
W1_hat_list <- purrr::map(res_list, "W1_hat")

# Save to file
S_string <- paste0("S", sum(groups))
for (i in 2:length(structure)) {
  S_string <- paste0(S_string, "-", length(structure[[i]]))
}
file <- paste0(
  S_string,
  "_T", T-h,
  "_M", M,
  "_dense"
)
saveRDS(results, file = paste("sim/sim_results/", file, ".rds", sep = ""))

saveRDS(error_list, file = paste("sim/sim_results/errorlist/", file, "_errorlist.rds", sep = ""))

saveRDS(W1_hat_list, file = paste("sim/sim_results/W1hat/", file, "_W1hat.rds", sep = ""))



# Inspect -----------------------------------------------


library(purrr)

## Line plot -------------------

MSE

MSE_ts <- transform_sim_MSE(MSE, F)

MSE_ts |> group_by(.model, h) |>
  summarise(mse = mean(MSE)) |>
  ggplot(aes(x = h, y = mse, color = .model)) +
  geom_line() +
  labs(x = "Horizon", y = "MSE") +
  theme_minimal()


MSE$mint_shr - MSE$mint_n

MSE_ts |>
  group_by(series, h) |>
  mutate(
    base_MSE = MSE[.model == "base"]
  ) |>
  ungroup() |>
  group_by(.model, h) |>
  summarise(
    MSE = mean(MSE),
    base_MSE = mean(base_MSE),
    pct_change = (MSE - base_MSE) / base_MSE * 100
  ) |>
  filter(h <=16) |>
  ggplot(aes(x = h, y = pct_change, color = .model)) +
  geom_line() +
  labs(x = "Horizon", y = "% improvements",
       title = "% relative improvements in MSE compared to Base") +
  theme_minimal()


## Box plot -------------------

# transform into df
error_df <- transform_error_list(error_list)

# box plot of 1-step-ahead error2
error_df %>%
  filter(h == 1) %>%  # filter for 1-step-ahead errors
  group_by(.model, id) %>%
  summarise(MSE = mean(e2)) %>%
  ggplot(aes(x = .model, y = MSE, color = .model)) +
    geom_boxplot() +
    theme_minimal()

# box plot of 1-step-ahead relative improvement
error_df %>%
  filter(h == 1) %>%
  group_by(.model, id) %>%
  summarise(MSE = mean(e2)) %>%
  # calculate relative improvement compared to base model
  group_by(id) %>%
  mutate(base_MSE = MSE[.model == "base"]) %>%
  ungroup() %>%
  mutate(pct_change = (MSE - base_MSE) / base_MSE * 100) %>%
  # remove outliers from mint
  filter(pct_change < 200) %>%
  # plot
  ggplot(aes(x = .model, y = pct_change, color = .model)) +
    geom_boxplot() +
    labs(x = "Model", y = "% relative improvements in MSE",
         title = "% relative improvements in MSE compared to base, 1-step-ahead forecasts") +
    theme_minimal()



