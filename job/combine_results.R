

path0 <- "job/S36-6-1/"
runid <- 2
struct <- "S36-6-1_"
metadata <- "T300_M100_"

filename <- paste0(struct, metadata)
path <- paste0(path0, "run", runid, "/")

# Initialise
template <- readRDS(paste0(path, filename, 1, ".rds"))
h <- dim(template$MSE[[1]])[1]

order_S <- rownames(template$S)

model_names <- c("base", "mint_shr", "mint_n", "mint_sample")
MSE_cum <- setNames(
  lapply(model_names, function(name) {
    matrix(0, h, length(order_S), dimnames = list(1:h, order_S))
  }),
  model_names
)
M <- 10000
W_shr_store <- numeric(M)
W_n_store   <- matrix(0, M, 2,
                      dimnames = list(NULL, c("lambda", "delta")))

# read files
for (i in 1:100) {
  fileread <- readRDS(paste0(path, filename, i, ".rds"))
  # index for 100 runs
  start <- (i - 1) * 100 + 1
  end <- i * 100
  W_shr_store[start:end] <- fileread$W_shr
  W_n_store[start:end, ] <- fileread$W_n
  for (model_name in model_names) {
    MSE_cum[[model_name]] <- MSE_cum[[model_name]] + fileread$MSE[[model_name]]
  }
}
MSE <- lapply(MSE_cum, function(x) x / 100)

# save results
template$MSE <- MSE
template$W_shr <- W_shr_store
template$W_n <- W_n_store

saveRDS(template, paste0(path0, metadata, "run", runid, ".rds"))


# ----------------------- errors_data

path_er <- paste0(path, "errors_data/")
er_list <- vector("list", 10000)
for (i in 1:100) {
  erread <- readRDS(paste0(path_er, filename, i, "_er.rds"))
  # index for 100 runs
  start <- (i - 1) * 100 + 1
  end <- i * 100
  er_list[start:end] <- erread
}

saveRDS(er_list, paste0(path0, metadata, "run", runid, "_er.rds"))
