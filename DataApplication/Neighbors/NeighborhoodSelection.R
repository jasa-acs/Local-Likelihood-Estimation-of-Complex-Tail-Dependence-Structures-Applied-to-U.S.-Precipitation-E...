### ========================================================================================================================================= ###
### Function to generated a homogeneous neighborhood (in terms of marginal tail distributions) around s0 based on the Hosking and Wallis test ###                              
### Daniela Castro-Camilo                                                                                                                     ###
### Email: daniela.castro.camilo@gmail.com                                                                                                    ###
### ========================================================================================================================================= ###

###################################
### Libraries and required data ###
###################################
# install.packages('rje')
library(rje)
library(fields)
load('DataApplication/Data/datamat5days1218.Rdata')
coord = read.table('DataApplication/Data/locations1218.txt', header = T); coord = coord[ , 2:3]
grid = read.table('DataApplication/Data/grid60.txt', header = T); grid = grid[, 2:3]

########################################################
### Function to generated a homogeneous neighborhood ###
########################################################
# data [matrix 2070x1218]: all the data (data) in original scale
# coord [matrix D0 x 2]: coordinates for all the stations
# grid [matrix 2235x3]: all grid points 
# s0: single grid point
# miles [logical]: should the distance be computed in miles? Default to false
# min.neigh [numeric]: minumum number of neighbors
# max.neigh [numeric]: maximum number of neighbors
# pr [numeric]: probability for the quantile-based treshold for each grid location
# alpha [numeric]: significance level for the tests
# dmax [numeric]: maximum distance (in km) where to look for neighbors
# rm.zeros [logical]: should zero precipitation be removed? Default to true
# which.test [vector]: which test? HW (1) or AD(2)

neighborhood_HT = function(data, coord, grid, s0, miles = FALSE, min.neigh = 5, max.neigh = 30, pr = 0.9, alpha = 0.05, dmax = 150, rm.zeros = TRUE, which.test = c(1,2)){
  source('DataApplication/Neighbors/HoTests.R')
  dist.mat = rdist.earth(coord, s0, miles = miles)
  distances = cbind(1:nrow(coord), dist.mat)
  distances = distances[order(distances[ , 2]), ]
  distances = distances[distances[, 2] <= dmax, ]
  if(!is.matrix(distances)){distances = matrix(distances, nrow = 1)}
  i = min(nrow(distances), max.neigh)
  if(i == 1){return(1)}
  if(nrow(distances) == 0){return(0)}
  
  repeat{
    cols = distances[1:i, 1]
    dists = distances[1:i, 2]
    eval.asses = assess.hom(data, cols, pr = pr, alpha = alpha, rm.zeros = rm.zeros, which.test = which.test)
    if(eval.asses$code == 1){
      n.neigh = length(cols)
      break
    }else{
      i = i - 1
    }
    if(i < 2) break
    
    if(exists("n.neigh")){
      if(n.neigh < min.neigh){cols = 0}else{cols = sort(distances[1:n.neigh, 1])}
    }else{
      cols = 0 
    }
  }
  return(cols)
}

# Takes time...
neigh = list()
for(i in 1:nrow(grid)){
  printPercentage(i, nrow(grid))
  s0 = grid[i, ]
  temp = tryCatch(neighborhood_HT(datamat, coord, grid, s0), error = function(e) e)
  if(!inherits(temp, "error"))
    neigh[[i]] = temp
}
save(neigh, file = 'Results/neigh_grid60km_dist150km.Rdata')
print('Done.')


