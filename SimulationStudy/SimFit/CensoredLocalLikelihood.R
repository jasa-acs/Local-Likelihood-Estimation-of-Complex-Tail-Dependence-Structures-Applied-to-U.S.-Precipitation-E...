### ==================================================================================================== ###
### Data generation, censored local log-likelihood, and model fitting functions for the simulation study ###
### Daniela Castro-Camilo                                                                                ###
### Email: daniela.castro.camilo@gmail.com                                                               ###
### ==================================================================================================== ###

#####################################
### Libraries and auxiliary codes ###
#####################################
# install.packages(c('mvtnorm', 'R.utils', 'fields', 'numDeriv', 'condMVNorm', 'Matrix', 'matrixcalc'))
library(mvtnorm); library(fields); library(numDeriv)
library(condMVNorm); library(Matrix); library(matrixcalc)
source('SimulationStudy/SimFit/Tools.R')
source('SimulationStudy/SimFit/logCopula.R')
source('SimulationStudy/SimFit/logLikelihood.R') # same as DataApplication/Fit/logLikelihood.R
source('SimulationStudy/SimFit/logPartial.R')

#######################
### Data simulation ###
#######################
# coord: locations
# N: number of samples
# neigh [numeric]: number of neighbors
# nu [numeric]: smoothing parameter
# b1, b2: define the different levels of non-stationarity used in Section 4 of the paper. See functions lambda.st() and rho() in SimulationStudy/SimFit/Tools.R for more details
# seed [numeric]: seed to specify with set.seed

data.sim = function(coord, N, neigh, nu, b1, b2, seed){
  nsim = nrow(coord)
  smooth = nu
  M <- Matern.st(coord, smothness = smooth, b1 = b1, b2 = b2)
  chol.M <- chol(M)
  lbda.st =  rep(0, nsim)
  for(i in 1:nsim)
    lbda.st[i] = lambda.st(as.numeric(coord[i, ]), smooth, b1, b2)
  
  sim = data.u = matrix(NA, N, nsim)
  set.seed(seed)
  for(i in 1:N){
    u = runif(1)
    e = (-1/lbda.st) * log(1 - u)
    sim[i, ] <- t(t(chol.M) %*% rnorm(nsim)) + e
  }
  for(i in 1:N){
    for(j in 1:nsim){
      w = sim[i,j]
      data.u[i,j] = F1(w, lbda.st[j])
    }
  }
  return(data.u)
}

##############################################
### Censored local log-likelihood function ###
##############################################
# theta [vector]: c(lambda, range), parameter vector
# data.u [matrix]: data matrix in uniform scale.
# coord [matrix D0 x 2]: coordinates for all the locations
# thres [numeric]: probability for the quantile-based treshold for each location
# nu [numeric]: fixed smoothing parameter
# censorL [logical]: should we compute censorized likelihood? Default to be true

model.likelihood = function(theta, data.u, coord, thres, nu, censorL){
  lbda = theta[1]; range = theta[2]; smooth = nu
  value1 = value2 = value3 = 0
  if(lbda <= 0 || range <= 0 ){
    return(1e09)
  } else{
    dist.mat = rdist(coord)# Distance matrix
    u.star = as.numeric(quantile(data.u[,1], probs = thres))
    # Classify the data# Classify the data into fully censored, non-censored, or partially censored
    N1 = dim(data.u)[1]; n1 = dim(data.u)[2]
    I = matrix(rep(1, N1 * n1), N1, n1)
    I[data.u <= u.star] = 0
    G1 = which(rowSums(I) == 0) # fully censored
    G2 = which(rowSums(I) == n1) # non-censored
    G3 = which(rowSums(I) > 0 & rowSums(I) < n1) # partially censored
    
    # Computing fully censored likelihood
    if(length(G1) > 0){
      if(censorL == TRUE)
        value1 <- -length(G1) * log.Cn(theta, u.star, dist.mat, nu)
      else
        value1 <- -log.lik(theta, data.u[G1, ], dist.mat, nu)
    }
    # Computing non-censored likelihood
    if(length(G2) > 0){
      value2 <- -log.lik(theta, data.u[G2, ], dist.mat, nu)
    }
    # Computing partially censored likelihood
    if(length(G3) > 0){
      if(censorL == TRUE)
        value3 <- -log.partialCn(theta, data.u[G3, ], u.star, dist.mat, nu)
      else
        value3 <- -log.lik(theta, data.u[G3, ], dist.mat, nu)
    }
  }
  # Negative log-likelihood
  return(sum(c(value1, value2, value3)))
}

##############################################################################################################
### Model fitting function: fits the model to the whole grid of size sqrt(grid.length) x sqrt(grid.length) ###
##############################################################################################################
# theta0 [vector 1x2]: c(lambda0, range0), initial values
# data.u [matrix]: data matrix in uniform scale
# coord [matrix]: coordinates for all the locations
# grid.length [matrix]: grid of size sqrt(grid.length) x sqrt(grid.length)
# thres [numeric]: probability for the quantile-based treshold for each location
# neigh [numeric]: number of neighbors
# nu [numeric]: fixed smoothing parameter
# censorL [logical]: should we compute censorized likelihood? Default to be true
# optimizer [character]: R function to optimize model.likelihood. Could be 'nlm' (recommended) or 'optim'

fit.model.likelihood = function(theta0, data.u, coord, grid.length, thres, neigh, nu, censorL, optimizer = 'nlm'){
  print('Fitting the model')
  # Prediction grid
  xgrid = ygrid = seq(1, 10, length.out = floor(sqrt(grid.length)))
  grid = expand.grid(xgrid, ygrid)
  # Optimization
  out_loc = matrix(NA, nrow(grid), 3)
  for(i in 1:nrow(grid)){
    t0 = grid[i, ]
    Nt0 = mi.kNN(t0, coord, neigh)
    data.t0 = data.u[ ,Nt0]
    coord.t0 = coord[Nt0, ]
    if(optimizer == 'nlm'){
      out <- nlm(model.likelihood, theta0, data.u = data.t0, coord = coord.t0, thres = thres, nu = nu, censorL = censorL)
      out_loc[i, 2:3] = out$estimate
    }
    if(optimizer == 'optim'){
      out <- optim(theta0, model.likelihood, data.u = data.t0, coord = coord.t0, thres = thres, nu = nu, censorL = censorL)
      out_loc[i, 2:3] = out$par 
    }
    out_loc[i, 1] = i
    out_loc = data.frame(out_loc)
    colnames(out_loc) = c('gridID', 'param1', 'param2')
  }
  return(out_loc)
}

