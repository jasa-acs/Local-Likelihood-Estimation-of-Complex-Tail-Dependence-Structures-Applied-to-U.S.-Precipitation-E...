### ============================================================================== ###
### Auxiliary functions to compute log copula function for partially censored data ###
### Daniela Castro-Camilo                                                          ###
### Email: daniela.castro.camilo@gmail.com                                         ###
### ============================================================================== ###

# theta = c(lambda, range)
# data.u: data matrix in uniform scale(n columns)
# u.star: treshold in u-scale
# dist.mat: distance matrix computed using rdist (package fields)
# nu: fixed smoothing parameter
dF1 <- function(z,theta,log=FALSE){
  lbda = theta[1]
  logdF1 <- log(lbda) -lbda*z + 0.5*lbda^2 + pnorm(z-lbda,log.p=TRUE)
  if(log){
    return(logdF1)
  } else{
    return(exp(logdF1))
  }
}

dFI <- function(z, I, dist.mat, theta, nu, log = FALSE){
  lbda = theta[1]
  range = theta[2]
  smooth = nu
  sigma <- Matern(dist.mat,range=range,smoothness=smooth)
  D <- ncol(dist.mat)
  k <- length(I)
  
  sigmaII <- sigma[I,I]
  C <- t(chol(sigmaII))
  Cinv <- solve(C)
  sigmaIIinv <- t(Cinv)%*%Cinv
  CinvzI <- Cinv%*%z[I]
  m1 <- t(CinvzI)%*%CinvzI
  m2 <- rep(1,k)%*%sigmaIIinv%*%z[I]
  m3 <- sum(sigmaIIinv)
  
  
  z.I <- rbind(0,z[-I]-sigma[-I,I]%*%sigmaIIinv%*%z[I])
  mu.I <- as.numeric((m2-lbda)/m3)*rbind(-1,1-sigma[-I,I]%*%sigmaIIinv%*%rep(1,k)) 
  sigma.I <- matrix(NA,D-k+1,D-k+1)
  sigma.I[1,1] <- m3^(-1)
  sigma.I[-1,1] <- (sigma[-I,I]%*%sigmaIIinv%*%rep(1,k)-1)/m3
  sigma.I[1,-1] <- t(sigma[-I,I]%*%sigmaIIinv%*%rep(1,k)-1)/m3
  sigma.I[-1,-1] <- sigma[-I,-I]-sigma[-I,I]%*%sigmaIIinv%*%sigma[I,-I]+(sigma[-I,I]%*%sigmaIIinv%*%rep(1,k)-1)%*%t(sigma[-I,I]%*%sigmaIIinv%*%rep(1,k)-1)/m3
  
  if (!isSymmetric(sigma.I))
    sigma.I = as.matrix(forceSymmetric(sigma.I))# warning("Problema numerico con la simetria de sigma.I Usamos forceSymmetric :/")
  
  logdFI <- log(lbda) -0.5*log(m3) -(k-1)*0.5*log(2*pi) -sum(log(diag(C))) -0.5*(m1 - (m2-lbda)^2/m3) + log(pmvnorm(upper=as.numeric(z.I),mean=as.numeric(mu.I),sigma=sigma.I))[1]
  if(log){
    return(logdFI)
  } else{
    return(exp(logdFI))
  }
}

dCI <- function(u, I, dist.mat, theta, nu, log = FALSE){
  z <- qF1(u,theta)
  logdCI <- dFI(z,I,dist.mat,theta,nu,log=TRUE) - sum(dF1(z[I],theta,log=TRUE))
  if(log){
    return(logdCI)
  } else{
    return(exp(logdCI))
  }
}
log.partialCn = function(theta, data.u, u.star, dist.mat, nu){
  value = NULL
  if( is.null(nrow(data.u)) ) data.u = matrix(data.u, nrow = 1, ncol = length(data.u)) # if data.u has only one row
  for(i in 1:nrow(data.u)){
    u = data.u[i, ]
    I = which(u > u.star)
    u[-I] = u.star
    value[i] = dCI(u, I, dist.mat, theta, nu, log = TRUE)
  }
  sum(value)
}
