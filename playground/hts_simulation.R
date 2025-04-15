k = 3
size = c(3, 5, 2)
rho = c(0.7, 0.7, 0.5)
delta = 0.39
epsilon = 0.99 - max(rho)
eidim = 5


ndim <- sum(size)
bigcor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
for (i in 1:k) {
  cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
  if (i == 1) {bigcor[1:size[1], 1:size[1]] <- cor}
  if (i != 1) {bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]),
                      (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor}
}
diag(bigcor) <- 1 - epsilon

eivect <- matrix(0, nrow = eidim, ncol = ndim)
for (i in 1:ndim) {
  ei <- runif(eidim, -1, 1)
  eivect[,i] <- sqrt(epsilon) * ei/sqrt(sum(ei^2))
}
bigE <- t(eivect) %*% eivect
cor.nz <- bigcor + bigE
cor.nz





