### ======================================================= ###
### Testing for Homogeneity of distributions among stations ###                              
### Daniela Castro-Camilo                                   ###
### Email: daniela.castro.camilo@gmail.com                  ###
### ======================================================= ###
# We test homogeneity of tail distributions among stations. Stations are ordered by distance. 
# If homogeneity is accepted for k stations, we proceed with station k+1 and so on until homogeneity is rejected.
# Specifically:
# Let Y_ij be the j-th observation at station i.
# Let X_ij = Y_ij/bar(Y), with bar(Y) a site specific index value (in this case, the median)
# If the observations are independent, we test H0: F_1 = ... = F_D = F, with F_i: distribution in the i-th station.
# F is a kappa(xi, alpha, k, h) distribution in the Hosking and Wallis test (h = 0 is GEV; h = 1 is GPD) and is left unspecified in the Anderson and Darling test.
# Reference: Viglione, A., Laio, F., Claps, P. (2007). "A comparison of homogeneity test for regional frequency analysis". Water resources research.

#####################################
### Libraries and auxiliary codes ###
#####################################
# install.packages('homtest')
library(homtest)
source("DataApplication/Neighbors/myHtests.R")

################################################
### Testing for Homogeneity of distributions ###
################################################
# data [matrix 2070x1218]: all the data (datamat) in original scale
# cols [vector]: column id to identify neighbors
# pr [numeric]: probability for the quantile-based treshold for each grid location
# alpha [numeric]: significance level for the tests
# rm.zeros [logical]: should zero precipitation be removed? Default to true
# DK.test [logical]: should we perform DK test?
# which.test [vector]: which test? HW (1) or AD(2)

assess.hom = function(data, cols, pr = 0.9, alpha = 0.05, rm.zeros = FALSE, DK.test = FALSE, which.test = c(1,2)){
  N.s0 = data[, cols]
  if(rm.zeros){N.s0[N.s0 == 0] = NA}
  z = NULL
  for(j in 1:length(cols)){
    y = N.s0[ , j]
    y = y[!is.na(y)]
    if(length(y) > 0){
      thres = quantile(y, pr)
      y = y[y > thres]
      temp = cbind(y, j)
      z = rbind(z, temp)
    }
  }
  x = z[ , 1]
  cod = z[ , 2]
  if(all(which.test == c(1,2))){
    test = tryCatch(myHW.tests(x, cod), error = function(e) e)
    if(!inherits(test, "error")){
      if(test[3] <= 0.23){
        value = test[1]
        if(value < 1){code = 1}else{code = 0} # 'acceptably homogeneous' if HW < 1
      } else{
        test2 = myADbootstrap.test(x, cod, Nsim = 500, alpha = alpha)
        if(DK.test) {test3 = myDK.test(x, cod)}
        if(test2[2] < (1 - alpha)){code = 1}else{code = 0}
        }
    } else{
      test2 = myADbootstrap.test(x, cod, Nsim = 500, alpha = alpha)
      if(test2[2] < (1 - alpha)){code = 1}else{code = 0}
      }
    } else{
      if(which.test == 1){
      test = myHW.tests(x, cod)
      value = test[1]
      if(value < 1){code = 1}else{code = 0}
    }
      if(which.test == 2){
      test2 = myADbootstrap.test(x, cod, Nsim = 500, alpha = alpha)
      if(test2[2] < (1 - alpha)){code = 1}else{code = 0}
    }
    }
  list(code = code, idcols = unique(cod))
}

