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



# groups <- c(2,2)
# groups <- c(4,4,4,4)
groups <- c(6,6,6,6,6,6)

T <- 316
h <- 16
Tsplit <- T - h

# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   list(c(1,2))
# )
# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   list(c(1,2), c(3,4)),
#   list(c(1,2))
# )
structure <- list(
  groups,
  as.list(seq(1,length(groups))),
  # list(c(1,2,3), c(4,5,6)),
  # list(c(1,2))
  list(1:6)
)

(S <- construct_S(
  structure = structure,
  sparse = FALSE,
  ascending = FALSE
))
order_S <- rownames(S)

# ranges for coefs in VAR
diag_range <- c(0.4, 0.7)
offdiag_range <- c(-0.4, 0.4)


A <- generate_block_diag(
  groups = groups,
  diag_range = diag_range,
  offdiag_range = offdiag_range,
  # message = message,
)$A

rho <- runif(length(groups), 0.6, 0.9)
Sigma <- generate_cor(
  groups = groups,
  rho = rho,
  delta = min(rho) * 0.5,
  # delta = 0.15,
  epsilon = (1-max(rho)) * 0.5,
  # epsilon = 0.15,
  eidim = length(groups)
)
# convert to cov using random sd
Sigma <- convert_cor2cov(Sigma)
# flip signs
V <- diag(x = sample(c(-1,1), size = sum(groups), replace = TRUE))
(Sigma <- V %*% Sigma %*% V)


# -------------------------------------------

params <- list(
  groups = groups,
  T = T,
  h = h,
  structure = structure,
  S = S,
  A = A,
  Sigma = Sigma
)

path <- "job/S36-6-1/run2/"
saveRDS(params, file = paste0(path, "params.rds"))

# -------------------------------------------



