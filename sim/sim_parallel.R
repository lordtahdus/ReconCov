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
groups <- c(4,4,4,4)
# groups <- c(6,6,6,6,6,6)

T <- 108
h <- 8
Tsplit <- T - h

# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   list(c(1,2))
# )
structure <- list(
  groups,
  as.list(seq(1,length(groups))),
  list(c(1,2,3,4))
)
# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   # list(c(1,2,3), c(4,5,6)),
#   # list(c(1,2))
#   list(1:6)
# )

(S <- construct_S(
  structure = structure,
  sparse = FALSE,
  ascending = FALSE
))
order_S <- rownames(S)

# ranges for coefs in VAR
diag_range <- c(0.4, 0.8)
offdiag_range <- c(-0.4, 0.4)

## VAR(1) block -------------------------
A <- generate_block_diag(
  groups = groups,
  diag_range = diag_range,
  offdiag_range = offdiag_range,
  stationary = FALSE,
)$A

plot_heatmap(A, TRUE)
A <- edit(A) ;colnames(A) <- NULL # edit manually

any((eigen(A)$values) < 0)


## Sigma ---------------------------------
rho <- runif(length(groups), 0.6, 0.8)
Sigma_raw <- generate_cor(
  groups = groups,
  rho = rho,
  delta = min(rho) * 0.5,
  # delta = 0.15,
  epsilon = (1-max(rho)) * 0.5,
  # epsilon = 0.15,
  eidim = length(groups)
)

plot_heatmap(Sigma_raw, TRUE)
Sigma <- edit(Sigma_raw) ; colnames(Sigma) <- NULL # edit manually
Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
any((eigen(Sigma)$values) < 1e-8)

# convert to cov using random sd
Sigma <- convert_cor2cov(Sigma)
# flip signs
V <- diag(x = sample(c(-1,1), size = sum(groups), replace = TRUE))
Sigma <- V %*% Sigma %*% V

plot_heatmap(Sigma %>% cov2cor(), TRUE)


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
  W_n <- novelist_cv(
    y,
    y_hat,
    S,
    window = round(Tsplit/2),
    message = FALSE
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

  # # # # # #
  # Return
  SSE <- list(
    base = ((actual - base_fc)^2),
    mint_shr = ((actual - recon_mint_shr)^2),
    mint_n = ((actual - recon_mint_n)^2),
    mint_sample = ((actual - recon_mint_sample)^2)
  )

  list(
    SSE = SSE,
    W_shr = W_shr$lambda,
    W_n = c(W_n$lambda, W_n$delta),
    W1_hat = y - y_hat
  )
}


# ---- PARALELL ----------------------------------
library(future.apply)
library(progressr)

handlers(global = TRUE) # Setup progress bar handler
handlers("txtprogressbar")  # or "progress" for a fancier bar

plan(multisession, workers = parallel::detectCores() - 2)

M <- 500

model_names <- c("base", "mint_shr", "mint_n", "mint_sample")
SSE_cum <- setNames(
  lapply(model_names, function(name) {
    matrix(0, h, length(order_S), dimnames = list(1:h, order_S))
  }),
  model_names
)

# PARALLEL
# res_list <- future_lapply(seq_len(M), function(i) run(), future.seed=TRUE)

# parallel with progress bar
with_progress({
  p <- progressor(along = 1:M)  # auto sets steps = length

  set.seed(1)
  res_list <- future_lapply(
    1:M, function(i) {
      p(message = sprintf("Sim %d", i))  # advances safely
    },
    future.seed=TRUE
  )
})

# check orders of results
res_list[[1]]$SSE |> names() == SSE_cum |> names()

# wrangle
SSE_cum <- Reduce(function(acc, res) Map(`+`, acc, res$SSE),
                      res_list, init = SSE_cum)
W_shr_store <- sapply(res_list, `[[`, "W_shr")
W_n_store <- t(sapply(res_list, `[[`, "W_n")) ; colnames(W_n_store) <- c("lambda", "delta")

MSE <- lapply(SSE_cum, function(mat) mat / M)

# plan(sequential) # Reset to sequential


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
  "_run3"
)
saveRDS(results, file = paste("sim/sim_results/", file, ".rds", sep = ""))

saveRDS(error_list, file = paste("sim/sim_results/", file, "_errorlist.rds", sep = ""))

saveRDS(W1_hat_list, file = paste("sim/sim_results/", file, "_W1hat.rds", sep = ""))

# Inspect --------------------


library(purrr)

MSE

MSE_ts <- transform_sim_MSE(MSE)

MSE_ts |> group_by(.model) |> index_by(h) |>
  summarise(mse = mean(MSE)) |>
  ggplot(aes(x = h, y = mse, color = .model)) +
  geom_line() +
  labs(x = "Horizon", y = "MSE") +
  theme_minimal()


MSE$mint_shr - MSE$mint_n


