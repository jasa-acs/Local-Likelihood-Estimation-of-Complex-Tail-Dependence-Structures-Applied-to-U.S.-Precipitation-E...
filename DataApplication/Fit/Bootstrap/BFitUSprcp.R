### ========================================================================== ###
### This is Fit/FitUsprcp.R adapted for the block-bootstrap samples            ###                              
### NOTE: this version runs on a single node in parallel across multiple cores ###
### Daniela Castro-Camilo                                                      ###
### Email: daniela.castro.camilo@gmail.com                                     ###
### ========================================================================== ###

#####################################
### Libraries and auxiliary codes ###
#####################################
library(parallel)
source('DataApplication/Fit/CensoredLocalLikelihood.R')

#######################################################
### Locations, grid, neighbors, and starting values ###
#######################################################
coord = read.table('DataApplication/Data/locations1218.txt', header = T) ; coord = coord[ , 2:3]
grid = read.table('DataApplication/Data/grid60.txt', header = T) 
load('DataApplication/Data/neigh_grid60km_dist150km.Rdata')
starting = read.table('DataApplication/Data/StartingValues.txt', header = T) 

#################################################################################
### Fit model to U.S. precipitation extremes for every block-Bootstrap sample ###
#################################################################################
# b [integer]: bootstrap sample
# rowID [numeric]: (As in Fit/FitUsprcp.R) row number in 'grid' that corresponds to the grid location where the model will be fitted 
# neigh [list]: (As in Fit/FitUsprcp.R) all neighbors for all grid points
# coord [matrix D0 x 2]: coordinates for all the stations
# grid [matrix 2235x3]: (As in Fit/FitUsprcp.R) all grid points where the model will be fitted
# starting [matrix 2235x2]: (As in Fit/FitUsprcp.R) c(lambda0, range0). Initial values
# thres [scalar]: (As in Fit/FitUsprcp.R) probability for the quantile-based treshold for each location. Default to be 0.8
# nu [scalar]: (As in Fit/FitUsprcp.R) fixed smoothing parameter. Default to be 0.5
# censorL [logical]: (As in Fit/FitUsprcp.R) should we compute censorized likelihood? Default to be true
# reltol [numeric]: (As in Fit/FitUsprcp.R) relative convergence tolerance. As in optim
# timeout [numeric]: (As in Fit/FitUsprcp.R) maximum number of seconds the function is allowed to run before being interrupted by the timeout. As in withTimeout
Bfit.USprcp <- function(rowID, b, neigh, coord, grid, starting, thres = 0.8, nu = 0.5, censorL = TRUE, reltol = 1e-3, timeout = 1800){
  load(paste0('DataApplication/Fit/Bootstrap/BootstrapSamples/datamatB', b,'.Rdata'))
  neighID = neigh[[rowID]]
  gridID = grid[rowID, ]
  theta0 = starting[rowID, ]
  if(all(!is.null(neighID) & neighID != 0)){
    out = fit.model.likelihood(gridID = gridID, data = new.datamat, theta0 = as.numeric(theta0), neigh = neighID, coord = coord, thres = thres, nu = nu, censorL = censorL, reltol = reltol, timeout = timeout)
  }else{
    out = "Empty neighborhood"
  }
  return(out)
}

#################################################################################################
### Example usage: fitting the model at grid locations 1 and 2 using block-Bootstrap sample 1 ###
#################################################################################################
print('Example usage: fitting the model at grid locations 1 and 2 using block-Bootstrap sample 1')
out = mclapply(1:2, Bfit.USprcp, b = 1, neigh = neigh, coord = coord, grid = grid, starting = starting, mc.cores = 2)
out[[1]]$par
out[[2]]$par
warning('Need to run DataApplication/Fit/Bootstrap/BlockBootstrap.R to create more bootstrap samples')

############################################
### Notes on parallelization of the code ###
############################################
# As mentioned above, this version of the code runs on a single node in parallel across multiple cores. But the number of grid points we can fit
# in parallel is limited to the the number of cores, mc.cores. 
# [Slow] To run the code for all grid points, fit.USprcp can be called sequentially using mclapply, e.g.,
# b = 1
# mc.cores = 20
# njobs = ceiling(nrow(grid)/mc.cores)
# out = list()
# for(i in 1:njobs){
#   resto = nrow(grid)%%mc.cores
#   tmp = if (i == njobs & resto > 0) resto else mc.cores
#   Rs = (mc.cores)*(i-1) + c(1:(tmp))
#   mc.cores = tmp
#   out[[i]] = mclapply(Rs, Bfit.USprcp, b = b, neigh = neigh, coord = coord, grid = grid, starting = starting, mc.cores = mc.cores)
# }
# [Faster] To run the code for all grid point when multiple nodes are available, it is advisable to submit each call of mclapply(Rs, fit.USprcp, ...)
# above to a different node.
