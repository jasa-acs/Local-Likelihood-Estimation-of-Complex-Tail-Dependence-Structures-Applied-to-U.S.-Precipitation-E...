### ========================================================================== ###
### Main function to fit the model to the simulated data                       ###                              
### NOTE: this version runs on a single node in parallel across multiple cores ###
### Daniela Castro-Camilo                                                      ###
### Email: daniela.castro.camilo@gmail.com                                     ###
### ========================================================================== ###

#####################################
### Libraries and auxiliary codes ###
#####################################
library(parallel)
source('SimulationStudy/SimFit/CensoredLocalLikelihood.R')

###################################
### Fit model to simulated data ###
###################################
fit.sim <- function(N, n, neigh, grid.length, starting, thres = 0.95, nu = 0.5, b1, b2, seed, censorL = TRUE){
  # Simulating locations
  nx = floor(sqrt(n))
  xx = seq(1, 10, length.out = nx)
  coord = expand.grid(xx,xx)
  # simulating data in uniform scale and fitting the model to each point in the sqrt(grid.length)xsqrt(grid.length) grid
  data.u = data.sim(coord, N, neigh, nu, b1, b2, seed)
  out = fit.model.likelihood(theta0 = starting, data.u = data.u, coord = coord, grid.length = grid.length, thres = thres, neigh = neigh, nu = nu, censorL = censorL)

  out
}
##############################################################################################
### Example usage: fitting the model to two simulated datasets (parallelizing over 'seed') ###
##############################################################################################
out = mclapply(12345:12346, fit.sim, N = 100, n = 100, neigh = 5, grid.length = 9, starting = c(0.5, 0.5),
               thres = 0.95, nu = 0.5, b1 = 5, b2 = 1.5, censorL = TRUE, mc.cores = 2)

out[[1]]
out[[2]]


############################################
### Notes on parallelization of the code ###
############################################
# As mentioned above, this version of the code runs on a single node in parallel across multiple cores. But the number of different 
# simulated datasets that we can fit in parallel is limited to the the number of cores, mc.cores. 
# [Slow] To run the code for Replic = 1000 datasets (as in the paper), fit.sim can be called sequentially using mclapply, e.g.,
# Replic = 1000
# mc.cores = detectCores()
# njobs = ceiling(Replic/mc.cores)
# out = list()
# for(i in 1:njobs){
#   resto = Replic%%mc.cores
#   tmp = if (i == njobs & resto > 0) resto else mc.cores
#   Rs = (mc.cores)*(i-1) + c(1:(tmp))
#   mc.cores = tmp
#   out[[i]] = mclapply(Rs, fit.sim, N = 100, n = 100, neigh = 5, grid.length = 9, starting = c(0.5, 0.5),
#                       thres = 0.95, nu = 0.5, b1 = 5, b2 = 1.5, censorL = TRUE, mc.cores = mc.cores)
# }
# [Faster] To run the code for Replic when multiple nodes are available, it is advisable to submit each call of mclapply(Rs, fit.sim, ...)
# above to a different node.


