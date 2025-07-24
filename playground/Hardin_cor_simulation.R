

##########################################################
#### R code for adding noise to correlation matrices  #### 
#### J. Hardin, S.R. Garcia, D. Golan                 ####
#### last updated 1/18/2013                           ####
##########################################################

#################
##   EXAMPLE   ##
#################


# adding noise to a 15x15 identity matrix

noise.iden <- noisecor(diag(15), epsilon = .5, eidim=2)


# adding noise to a 2 block constant correlation matrix
rho.blck = c(.9,.2)
eps.blck = 0.99 - max(rho.blck)
samp.blck <- simcor(k = 2, size = c(10, 5), rho = c(0.9, 0.2), delta = 0.39, epsilon = eps.blck, eidim = 2) 


# adding noise to a 2 block correlation matrix with an AR(1) strucutre
rho.Top = c(.7,.9)
eps.Top = (1-max(rho.Top))/(1+max(rho.Top)) - .01
samp.Top <- simcorTop(k=2, size=c(10,5), rho=c(.7,.9), epsilon = eps.Top, eidim=2)


# adding noise to a 2 block correlation matrix with a Hub structure
rho.Hub = c(.9,.7)
tau.Hub = c( (.9 - .7) / (10-2), (.7-.6)/(5-2))
eps.Hub = min(1-(rho.Hub) - 0.75*(tau.Hub) ) - .01
samp.Hub <- simcor.H(k=2, size=c(10,5),	rho=rbind(c(.9,.7), c(.7,.6)), power=1,
                     epsilon=eps.Hub, eidim=2)



########################
## The Basic Function ##
########################

# this adds noise to any user specified correlation matrix

noisecor <- function(cormat, epsilon = .01, eidim=2){
  
  ndim=dim(cormat)[1]
  diag(cormat) <- 1 - epsilon
  
  ### adding noise to the correlation matrix
  
  eivect <- c( )
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
  }
  
  
  bigE <- t(eivect) %*% eivect
  cor.nz <- cormat + bigE
  cor.nz
  
}


#######################################
##  Simulating different correlation ##
##  structures as in Hardin et al.   ##
#######################################


################################
##  Simulating the Constant Correlation
################################

# this function simulates a block correlation matrix
# The size and base correlation for each block is user specified
# There is an additional delta parameter for the off diagonal correlations


# k is the number of groups
# size is a vector of length k specifying the size of each group 
# rho is a vector of length k specifying base correlation values
# epsilon <- 0.99 - max(rho)
# eidim is the space from which the noise is generated, the smaller the more noise
# delta is the correlation of the off diagonal blocks

simcor = function (k = 6, size = c(10, 5, 8, 2, 15, 50), rho = c(0.7, 
                                                                 0.7, 0.5, 0.9, 0.85, 0.4), delta = 0.39, epsilon = 0.99 - 
                     max(rho), eidim = 2) 
{
  ndim <- sum(size)
  bigcor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
  for (i in 1:k) {
    cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
    if (i == 1) {bigcor[1:size[1], 1:size[1]] <- cor}
    if (i != 1) {bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]), 
                        (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  
  eivect <- c()
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  cor.nz
}


################################
##  Simulating the Toeplitz Matrix
################################

# this function calculates an AR(1) Toeplitz matrix
# with a block diagonal structure.  The size and base correlation for each
# block is user specified


# k is the number of groups
# size is a vector of length k specifying the size of each group 
# rho is a vector of length k specifying base correlation values
# epsilon <- (1-max(rho))/(1+max(rho)) - .01 
# eidim is the space from which the noise is generated, the smaller the more noise

simcorTop <- function(k=6,size=c(10,5,8,2,15,50),rho=c(.7,.7,.5,.9,.85,.4),epsilon= .01, eidim=2) {
  
  
  ndim <- sum(size)# dim of correlation matrix
  bigcor<- matrix(rep(0, ndim*ndim), ncol=ndim)
  
  ### generating the basic correlation matrix
  
  
  for (i in 1:k){
    
    top <- c(1,rho[i]^(seq(1:(size[i]-1))))
    cor <- toeplitz(top)
    
    if (i==1){bigcor[1:size[1], 1:size[1]] <- cor}
    if (i!=1){bigcor[(sum(size[1:(i-1)]) + 1):sum(size[1:i]),
                     (sum(size[1:(i-1)]) + 1):sum(size[1:i])] <- cor}
  }
  
  diag(bigcor) <- 1 - epsilon
  
  
  ### adding noise to the correlation matrix
  
  eivect <- c( )
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
  }
  
  
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  cor.nz
  
}

################################
##  Simulating the Hub Matrix (entries filled in using Toeplitz structure)
################################

# this function calculates a Toeplitz matrix with values descending
# from a user specified maximum to minimum.  The matrix has a 
# block diagonal structure.  The size and base correlation for each
# block is user specified.


# k is the number of groups
# size is a vector of length k specifying the size of each group 
# rho is a vector of length k specifying base correlation values
# epsilon <- (1-min(rho) - 0.75*min(tau) ) - .01
# tau_k = (max(rho_k) - min(rho_k) )/ (size_k -2) 
# eidim is the space from which the noise is generated, the smaller the more noise
# power = 2 makes the correlations stay high
# power = 0.5 makes the correlations descent rapidly


simcor.H <- function(k=6, size=c(10,5,8,7,15,50), 
                     rho=rbind(c(.9,.7), c(.7,.7), c(.7,.2), c(.5,.3), c(.9,.85), c(.3,.2)), power=1,
                     epsilon=.08, eidim=2){
  
  
  ndim <- sum(size)# dim of correlation matrix
  bigcor<- matrix(rep(0, ndim*ndim), ncol=ndim)
  
  ### generating the basic correlation matrix
  
  
  for (i in 1:(k) ){
    
    cor <- toeplitz(rho.func(rho[i,1],rho[i,2],power,size[i]) )
    
    if (i==1){bigcor[1:size[1], 1:size[1]] <- cor}
    if (i!=1){bigcor[(sum(size[1:(i-1)]) + 1):sum(size[1:i]),
                     (sum(size[1:(i-1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  
  
  ### adding noise to the correlation matrix
  
  eivect <- c( )
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
  }
  
  
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  cor.nz
}


# rho.func is needed for filling in the rest of the structure of the Hub
# correlation matrix
# r.max is the maximum user specified correlation
# r.min is the minimum user specified correlation
# power is the power at which the correlations descend
# p is the size of the correlation block

rho.func <- function(r.max, r.min, power,p){
  rhovec <-c()
  
  rhovec[1] <- 1
  for(i in 2:p){
    rhovec[i] <- r.max - ((i-2)/(p-2))^power*(r.max-r.min)
  }
  rhovec}



