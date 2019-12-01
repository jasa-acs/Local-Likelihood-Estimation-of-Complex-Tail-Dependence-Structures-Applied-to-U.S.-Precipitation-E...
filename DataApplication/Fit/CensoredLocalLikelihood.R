### ========================================================================================================= ###
### Censored local log-likelihood and model fitting functions for the U.S. precipitation extremes application ###
### Daniela Castro-Camilo                                                                                     ###
### Email: daniela.castro.camilo@gmail.com                                                                    ###
### ========================================================================================================= ###

#####################################
### Libraries and auxiliary codes ###
#####################################
# install.packages(c('mvtnorm', 'R.utils', 'fields', 'numDeriv', 'condMVNorm', 'Matrix', 'matrixcalc'))
library(mvtnorm); library(R.utils); library(fields); library(numDeriv)
library(condMVNorm); library(Matrix); library(matrixcalc)
source('DataApplication/Fit/Tools.R')
source('DataApplication/Fit/logCopula.R')
source('DataApplication/Fit/logLikelihood.R')
source('DataApplication/Fit/logPartial.R')

##############################################
### Censored local log-likelihood function ###
##############################################
# theta [vector]: c(lambda, range), parameter vector
# data.u [matrix 2070 x D0]: data matrix in uniform scale. D0 is the number of neighbors
# coord [matrix D0 x 2]: coordinates for all the stations
# thres [numeric]: probability for the quantile-based treshold for each location
# nu [numeric]: fixed smoothing parameter
# censorL [logical]: should we compute censorized likelihood? Default to be true

model.likelihood = function(theta, data.u, coord, thres, nu, censorL){
  lbda = theta[1]; range = theta[2]; smooth = nu
  value1 = value2 = value3 = 0
  if(lbda <= 0 || range <= 0 ){
    return(1e09)
  }else{
    dist.mat = rdist.earth(coord, miles = F)
    u.star = apply(data.u, 2, quantile, probs = thres, na.rm = T)
    # Classify the data into fully censored, non-censored, or partially censored
    N1 = dim(data.u)[1]; n1 = dim(data.u)[2]
    I = matrix(rep(NA, N1 * n1), N1, n1)
    for(i in 1:N1){
      cond0 = data.u[i, ] <= u.star
      cond1 = data.u[i, ] > u.star
      I[i, which(cond0)] = 0
      I[i, which(cond1)] = 1
    }
    G1 = which(rowSums(I, na.rm = T) == 0) # fully censored
    G2 = which(rowSums(I, na.rm = T) == n1) # non-censored
    G3 = which(rowSums(I, na.rm = T) > 0 & rowSums(I, na.rm = T) < n1) # partially censored
    
    # Computing fully censored likelihood
    if(length(G1) > 0){
      if(censorL == TRUE){
        if(all(!is.na(data.u[G1, ]))){
          value1 = -length(G1) * log.Cn(theta, u.star, dist.mat, nu)
        }else{
          ## Handling rows with NAs
          row.na = NULL
          f.na = function(x) any(is.na(x))
          data.u.g1 = data.u[G1, ]
          if(is.null(nrow(data.u.g1))) data.u.g1 = matrix(data.u.g1, nrow = 1, ncol = length(data.u.g1))
          for(k in 1:nrow(data.u.g1)){
            x = data.u.g1[k, ]
            if(f.na(x)) row.na = c(row.na, k)
          }
          value1.j = NULL
          for(j in row.na){
            y.j = data.u.g1[j, ]
            id.na = which(is.na(y.j))
            dist.mat.j = dist.mat[-id.na, -id.na]
            u.star.j = u.star[-id.na]
            value1.j = c(value1.j, -log.Cn(theta, u.star.j, dist.mat.j, nu))
          }
          value11 <- sum(value1.j)
          ## Handling rows without NAs
          row.nan = (1:length(G1))[-row.na]
          value12 = 0
          if(length(row.nan) > 0) value12 <- -length(row.nan) * log.Cn(theta, u.star, dist.mat, nu)
          value1 = value11 + value12
        }
      }
      else
        value1 = -log.lik(theta, data.u[G1, ], dist.mat, nu)
    }
    # Computing non-censored likelihood
    if(length(G2) > 0){
      value2 = -log.lik(theta, data.u[G2, ], dist.mat, nu)
    }
    # Computing partially censored likelihood
    if(length(G3) > 0){
      if(censorL == TRUE)
        value3 = -log.partialCn(theta, data.u[G3, ], u.star, dist.mat, nu)
      else
        value3 = -log.lik(theta, data.u[G3, ], dist.mat, nu)
    }
    # Negative log-likelihood
    return(sum(c(value1, value2, value3))) 
  }
}

##############################
### Model fitting function ###
##############################
# gridID [vector 1x3]: c(id, long, lat). Location of the grid where we want to fit our model
# data [matrix 2070x1218]: all the data (datamat) in original scale
# theta0 [vector 1x2]: c(lambda0, range0), initial values
# neigh [vector]: id of the neighbors of gridID. It corresponds to column numbers in 'data'
# coord [matrix D0 x 2]: as in model.likelihood
# thres [numeric]: as in model.likelihood
# nu [numeric]: as in model.likelihood
# censorL [logical]: as in model.likelihood
# reltol [numeric]: relative convergence tolerance. As in optim
# timeout [numeric]: maximum number of seconds the function is allowed to run before being interrupted by the timeout. As in withTimeout

fit.model.likelihood = function(gridID, data, theta0, neigh, coord, thres, nu, censorL = TRUE, reltol, timeout){
  coord = coord[neigh, ]
  datamat.s0 = data[ , neigh]
  datamat.s0[datamat.s0 == 0] = NA # we remove zeros
  # If all the row is NA
  rm.idx = NULL
  f.na = function(x) all(is.na(x))
  for(i in 1:nrow(datamat.s0)){
    if(f.na(datamat.s0[i, ])) rm.idx = c(rm.idx, i)
  }
  if(!is.null(rm.idx)) datamat.s0 = datamat.s0[-rm.idx, ]
  # If all the row is NA, except one
  rm(rm.idx); rm.idx = NULL
  for(i in 1:nrow(datamat.s0)){
    x = datamat.s0[i, ]
    if(sum(!is.na(x)) == 1){
      rm.idx = c(rm.idx, i)
    }
  }
  if(!is.null(rm.idx)) datamat.s0 = datamat.s0[-rm.idx, ]
  
  data.u = apply(datamat.s0, 2, u.np) # transform to uniform margins
  out = list()
  dummyFUN = function(theta0, data.u, coord, thres, nu, censorL, trace){
    temp = tryCatch(optim(theta0, model.likelihood,  data.u = data.u, coord = coord, thres = thres, nu = nu,
                          censorL = censorL, control = list(reltol = reltol)), error = function(e) e)
    if(!inherits(temp, "error") & is.numeric(temp$par)){
      out = temp; out$start = theta0; out$gridID = gridID
    }
    out
  }
  results = tryCatch(withTimeout({dummyFUN(theta0, data.u, coord, thres, nu, censorL, trace)}, timeout = timeout, onTimeout = "warning"), 
                     error = function(e) e)
  if(!inherits(results, "error")){
    output = results
  }else{
    output = gridID
  }
  return(output)
}

