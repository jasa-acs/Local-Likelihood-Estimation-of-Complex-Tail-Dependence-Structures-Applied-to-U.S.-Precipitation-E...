### =========================================================================== ###
### Auxiliary function to compute log likelihood function for non-censored data ###
### Daniela Castro-Camilo                                                       ###
### Email: daniela.castro.camilo@gmail.com                                      ###
### =========================================================================== ###

# theta = c(lambda, range)
# data.u: data matrix in uniform scale(n columns)
# dist.mat: distance matrix computed using rdist (package fields)
# nu: fixed smoothing parameter
log.lik = function(theta, data.u, dist.mat, nu){
  lbda = theta[1]; range = theta[2]; smooth = nu
  # lbda = theta[1]; range = theta[2]; smooth = nu
  if( is.null(nrow(data.u)) ) data.u = matrix(data.u, nrow = 1, ncol = length(data.u)) # if data.u has only one row
  z = apply(data.u, 1:2, F1inv, lbda = lbda)
  # z = ToolsPkg::F1inv_apply(data.u, lbda, 1e-09, 1e03)
  # Matern covariance
  sigma = Matern(dist.mat, range = range, smoothness = smooth)
  n = ncol(z)
  C = t(chol(sigma)) 
  C.inv = solve(C)
  sigma.inv = t(C.inv) %*% C.inv
  
  log.f1 = function(zij, lbda){
    if(is.na(zij)) return(NA)
    log(lbda) + lbda^2/2 - lbda * zij + pnorm(zij - lbda, log.p = T)
  }
  
  log.fn = function(zi, lbda, sigma, sigma.inv, C.inv){
    ## To remove observations with NAs
    if(any(is.na(zi))){
      id.na = which(is.na(zi))
      zi = zi[!is.na(zi)]
      n = length(zi)
      sigma = sigma[-id.na, -id.na]
      sigma.inv = sigma.inv[-id.na, -id.na]
      C.inv = C.inv[-id.na, -id.na]
    }
    ##
    m1.1 = zi%*%t(C.inv)
    m1 = t(m1.1)%*%m1.1
    m1 = t(zi) %*% sigma.inv %*% zi
    m2 = t(rep(1, n)) %*% sigma.inv %*% zi
    m3 = sum(sigma.inv)
    m1.star = (m2 - lbda)/sqrt(m3)
    out = log(lbda)  - (n - 1)/2 * log(2 * pi) - sum(log(diag(C))) - (1/2) * log(m3) + ((m1.star) ^ 2 - m1)/2 + pnorm(m1.star, log.p = T)
    return(out)
  }
  
  lfn = apply(z, 1, log.fn, lbda, sigma, sigma.inv, C.inv)
  lf1 = apply(z, c(1,2), log.f1, lbda)
  value <- sum(lfn) - sum(lf1, na.rm = T)
  value
}