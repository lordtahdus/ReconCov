library(MASS)
# library(matrixcalc)
library(Matrix)
library(tidyr)

library(fabletools)
library(fable)
library(feasts)
library(tsibble)
library(dplyr)

library(devtools)
load_all()


# Parameters -----------------------------------
M <- 100
path <- "job/S16-4-1/run2/"

params <- readRDS(paste0(path, "params.rds"))

groups <- params$groups
A <- params$A
Sigma <- params$Sigma
T <- params$T
h <- params$h
structure <- params$structure
S <- params$S

Tsplit <- T - h
order_S <- rownames(S)



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
    W_n = c(W_n$lambda, W_n$delta)
  )
}




model_names <- c("base", "mint_shr", "mint_n", "mint_sample")
SSE_cum <- setNames(
  lapply(model_names, function(name) {
    matrix(0, h, length(order_S), dimnames = list(1:h, order_S))
  }),
  model_names
)

W_shr_store <- numeric(M)
W_n_store   <- matrix(0, M, 2,
                      dimnames = list(NULL, c("lambda", "delta")))
error_list <- vector("list", M)

for(i in seq_len(M)) {
  cat("Iteration ", i, "\n")

  res <- run(A, Sigma)
  # add up the SSE matrices
  # store the error list
  # SSE_cum2 <- Map(`+`, SSE_cum2, res$SSE)
  error_list[[i]] <- res$SSE
  # store the parameters
  W_shr_store[i]         <- res$W_shr
  W_n_store[i, ]         <- res$W_n
}


SSE_cum <- setNames(
  lapply(model_names, function(name)
    Reduce(`+`, lapply(error_list, `[[`, name))
  ),
  model_names
)
# calculate MSE from error list
MSE <- lapply(SSE_cum, function(mat) mat / M)


# Save --------------------------------

# Combine all your outputs into a named list
sim_results <- list(
  groups = groups,
  S      = S,
  A      = A,
  Sigma  = Sigma,
  MSE    = MSE,    # or averaged MSE
  W_shr  = W_shr_store,       # can be a list of matrices or one matrix
  W_n    = W_n_store          # same
)

# Save to file
args <- commandArgs(trailingOnly = TRUE)
# start <- as.integer(args[1])
# end <- as.integer(args[2])
index <- as.integer(args[1])

S_string <- paste0("S", sum(groups))
for (i in 2:length(structure)) {
  S_string <- paste0(S_string, "-", length(structure[[i]]))
}

file <- paste0(
  S_string,
  "_T", T-h,
  "_M", M,
  # start, "_", end
  "_", index
)
saveRDS(sim_results, file = paste(path, file, ".rds", sep = ""))
saveRDS(error_list, file = paste(path,"errors_data/",file, "_er.rds", sep = ""))







