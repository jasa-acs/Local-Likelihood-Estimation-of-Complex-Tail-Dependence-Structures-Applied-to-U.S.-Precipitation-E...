### ========================================================================== ###
### Main function to fit the model to U.S. precipitation extremes              ###                              
### NOTE: this version runs on a single node in parallel across multiple cores ###
### Daniela Castro-Camilo                                                      ###
### Email: daniela.castro.camilo@gmail.com                                     ###
### ========================================================================== ###

#####################################
### Libraries and auxiliary codes ###
#####################################
library(parallel)
source('DataApplication/Fit/CensoredLocalLikelihood.R')

#############################################################
### Data, locations, grid, neighbors, and starting values ###
#############################################################
load('DataApplication/Data/datamat5days1218.Rdata')
coord = read.table('DataApplication/Data/locations1218.txt', header = T) ; coord = coord[ , 2:3]
grid = read.table('DataApplication/Data/grid60.txt', header = T) 
load('DataApplication/Data/neigh_grid60km_dist150km.Rdata')
starting = read.table('DataApplication/Data/StartingValues.txt', header = T) 
# NOTE: 
# Starting values were obtained from a first run of the model, and are used to speed up the code. 
# As an alternative, suitable random starting values can be given, e.g., as
# starting = matrix(c(runif(nrow(grid), 2, 8), runif(nrow(grid), 10, 50)), nrow = nrow(grid))

################################################
### Fit model to U.S. precipitation extremes ###
################################################
# rowID [numeric]: row number in 'grid' that corresponds to the grid location where the model will be fitted 
# data [matrix 2070x1218]: all the data (datamat) in original scale
# neigh [list]: all neighbors for all grid points.
# coord [matrix D0x2]: coordinates for all locations.
# grid [matrix 2235x3]: all grid points where the model will be fitted
# starting [matrix 2235x2]: c(lambda0, range0). Initial values
# thres [scalar]: probability for the quantile-based treshold for each location. Default to be 0.8
# nu [scalar]: fixed smoothing parameter. Default to be 0.5
# censorL [logical]: should we compute censorized likelihood? Default to be true
# reltol [numeric]: relative convergence tolerance. As in optim
# timeout [numeric]: maximum number of seconds the function is allowed to run before being interrupted by the timeout. As in withTimeout
fit.USprcp <- function(rowID, data, neigh, coord, grid, starting, thres = 0.8, nu = 0.5, censorL = TRUE, reltol = 1e-3, timeout = 1800){
  neighID = neigh[[rowID]]
  gridID = grid[rowID, ]
  theta0 = starting[rowID, ]
  if(all(!is.null(neighID) & neighID != 0)){
    out = fit.model.likelihood(gridID = gridID, data = data, theta0 = as.numeric(theta0), neigh = neighID, coord = coord, thres = thres, nu = nu, censorL = censorL, reltol = reltol, timeout = timeout)
  }else{
    out = "Empty neighborhood"
  }
  return(out)
}

##################################################################
### Example usage: fitting the model at grid locations 1 and 2 ###
##################################################################
print('Example usage: fitting the model at grid locations 1 and 2')
out = mclapply(1:2, fit.USprcp, data = datamat, neigh = neigh, coord = coord, grid = grid, starting = starting, mc.cores = 2)
out[[1]]$par
out[[2]]$par

############################################
### Notes on parallelization of the code ###
############################################
# As mentioned above, this version of the code runs on a single node in parallel across multiple cores. But the number of grid points we can fit
# in parallel is limited to the the number of cores, mc.cores. 
# [Slow] To run the code for all grid points, fit.USprcp can be called sequentially using mclapply, e.g.,
# mc.cores = detectCores()
# njobs = ceiling(nrow(grid)/mc.cores)
# out = list()
# for(i in 1:njobs){
#   resto = nrow(grid)%%mc.cores
#   tmp = if (i == njobs & resto > 0) resto else mc.cores
#   Rs = (mc.cores)*(i-1) + c(1:(tmp))
#   mc.cores = tmp
#   out[[i]] = mclapply(Rs, fit.USprcp, data = datamat, neigh = neigh, coord = coord, grid = grid, starting = starting, mc.cores = mc.cores)
# }
# [Faster] To run the code for all grid point when multiple nodes are available, it is advisable to submit each call of mclapply(Rs, fit.USprcp, ...)
# above to a different node.
