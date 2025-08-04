library(Matrix)
library(dplyr)

library(devtools)
load_all()

resid <- y - y_hat

# choosing number of factors K
eigenvals <- eigen(cov(resid))$values 
eigenvals |> plot()
sum(eigenvals > 1)

# Bai & Ng (2002) Information Criterion
# ...

library(nFactors)
nScree(x = eigenvals)$Components
nScree(x = eigenvals) |> plotnScree()

# factors
K <- 2
pca <- prcomp(resid)
B <- pca$rotation[,1:K] # loadings
F <- pca$x[,1:K] # factor scores

remainder <- resid - F %*% t(B) # residuals

all.equal(
  B %*% diag(pca$sdev[1:K]^2) %*% t(B) + cov(remainder),
  cov(resid),
  tolerance = 0
)

factor_analysis <- factanal(resid, factors = K, rotation = "none")
factor_analysis$loadings

# Cov ests
W_n <- novelist_cv(
  y, y_hat,
  S, 
  window_size = 60, 
  h = 1,
)
W_n$lambda
W_n$delta

eigen(cov(W_n$cov))$values

W_shr <- shrinkage_est(
  resid
)
W_shr$lambda
