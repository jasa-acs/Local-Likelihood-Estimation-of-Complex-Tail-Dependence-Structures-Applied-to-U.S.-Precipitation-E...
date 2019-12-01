### ========================================================================= ###
### Auxiliary function to compute log copula function for fully censored data ###
### Daniela Castro-Camilo                                                     ###
### Email: daniela.castro.camilo@gmail.com                                    ###
### ========================================================================= ###

# theta = c(lambda, range)
# u.star: threshold (vector) in uniform scale
# dist.mat: distance matrix computed using rdist (package fields)
# nu: fixed smoothing parameter
log.Cn = function(theta, u.star, dist.mat, nu){ 
  # Warning: works only for "sigma" such that diag(sigma) = 1
  lbda = theta[1]; range = theta[2]; smooth = nu
  z = sapply(u.star, F1inv, lbda = lbda) # = z.star
  sigma = Matern(dist.mat, range = range, smoothness = smooth)
  z = rep(z, ncol(sigma))
  n = ncol(sigma)
  mean = rep(0, n)
  probs = NULL
  condVars = lapply(1:n, mi_condMVN, mean = mean, sigma = sigma)
  unos = rep(1, (n - 1))
  sigma0 = matrix(NA, nrow = n, ncol = n)
  sigma0[n, n] = 1
  for(i in 1 : n){
    tmp = condVars[[i]]
    condVar = tmp$condVar
    C = tmp$C
    D = (unos - C); Dt = t(D)
    sigma0[1:(n - 1), 1:(n - 1)] = condVar + D %*% Dt
    sigma0[1:(n - 1), n] = -D
    sigma0[n, 1:(n - 1)] = -Dt
    upper = c(z[-i] - z[i]*unos + D * lbda, z[i] - lbda)
    set.seed(302)
    probs[i] = exp(lbda ^ 2/2 - lbda * z[i]) *pmvnorm(lower = -Inf, upper = upper, sigma = sigma0, algorithm = GenzBretz())[1]
  }
  set.seed(302)
  value = as.numeric(pmvnorm(lower = -Inf, upper = z, sigma = sigma)[1] - sum(probs))
  log(value)
}
