### ========================================================== ###
### Auxiliary functions to fit the model to the simulated data ###
### Daniela Castro-Camilo                                      ###
### Email: daniela.castro.camilo@gmail.com                     ###
### ========================================================== ###

#################
### Functions ###
#################
# F.rank: transform to uniform scale
# F1: computed the univarite cdf for the stationary factor copula model.
# dummy: auxiliar function to compute the inverse of F1.
# F1inv: the inverse of F1.
# mi.kNN: k nearest neighboors
# mi_condMVN : Computes conditional mean and covariance to be used in logCopula.R
# lambda.st:  non-stationary rate
# rho: non-stationary range
# Mnu, Kst, Matern.st: to compute non-stationary Matern covariance function 

F.rank = function(x){
  N = length(x)
  (rank(x)-0.5)/N
}

F1 = function(w, lbda){
  pnorm(w) - pnorm(w - lbda) * exp(lbda^2/2 - lbda * w)
}

dummy = function(w, u, lbda){
  return(F1(w, lbda) - u)
}

F1inv = function(u, lbda, tol = 1e-09, niter = 1e03){
  x0 = qnorm(u, 0, 1) - 0.5
  x1 = qnorm(u, 0, 1)
  for(i in 1: niter){
    temp0 = dummy(x0, u, lbda)
    temp1 = dummy(x1, u, lbda)
    x2 = x1 - temp1 * (x1 - x0)/(temp1 - temp0)
    temp2 = dummy(x2, u, lbda)
    aux2 = abs(temp2)
    aux21 = abs(x2-x1)
    if(aux2 < tol) return(x2)
    x0 = x1
    x1 = x2
  }
  return(x2)
}

qF1 = function(u, theta){
  lbda = theta[1]
  return(sapply(u, F1inv, lbda = lbda))
}


mi.kNN = function(x1, x2, k) { # x1 is a vector; x2 is a nx2 matix
  n = dim(x2)[1]
  # k = k + 1 # shift to avoid outputing x1
  temp2 = NULL
  for(i in 1:n){
    temp1 = x2[i,]
    temp2[i] = sqrt(sum((x1 - temp1) ^ 2))
  }
  temp3 = cbind(1:n,temp2)
  temp3 = temp3[order(temp3[,2]),]
  idx = temp3[1:k,1]
  d.idx = as.numeric(temp3[1:k,2])
  knn = x2[idx,]
  idx
}

mi_condMVN = function(mean, sigma, given.ind){
  dependent.ind = (1:length(mean))[-given.ind]
  B = sigma[dependent.ind, dependent.ind]
  C = sigma[dependent.ind, given.ind, drop = FALSE]
  D = sigma[given.ind, given.ind]
  CDinv = C %*% solve(D)
  Dinv = solve(D)
  cVar = B - CDinv %*% t(C)
  list(C = C, condVar = cVar)
}


lambda.st = function(s, nu, b1, b2){
  theta.st = function(s)
    (1/9) * (0.6 * s[2] + 10.2)
  
  t = c(s[1] + 1, s[2])
  r = Kst(s, t, nu, b1, b2)
  qnorm(theta.st(s)/2) * (2*(1 - r))^(-1/2)
}

rho = function(s, b1, b2)
  1 - pnorm(s[1], b1, b2) + 0.5

Mnu = function(arg, nu)
  abs(arg) ^ nu * besselK(arg, nu)

Kst = function(s, t, smothness, b1, b2){
  if((s[1] == t[1]) & (s[2] == t[2]))
    return(1)
  nu = smothness
  tmp = rep(NA, 3)
  d = as.numeric(dist(rbind(s, t)))
  if (any(d < 0)) 
    stop("distance argument must be nonnegative")
  d[d == 0] <- 1e-10
  d = d/(2*sqrt(nu)) ###### Para hacer coincidir con Matern() de package fields
  con = (2^(nu - 2)) * gamma(nu)
  tmp[1] = 1/con
  den = rho(s, b1, b2) ^ 2 + rho(t, b1, b2) ^ 2
  tmp[2] = (rho(s, b1, b2) * rho(t, b1, b2))/den
  arg = (2 * sqrt(2 * nu) * d)/sqrt(den)
  tmp[3] = Mnu(arg, nu)
  return(prod(tmp))
}

Matern.st = function(loc, smothness, b1, b2){
  nsim = nrow(loc)
  pares = t(combn(nsim, 2))
  matern.mat = matrix(NA, nsim, nsim)
  diag(matern.mat) = 1
  for(i in 1:nrow(pares)){
    s1 = pares[i,1]
    s2 = pares[i,2]
    s = as.numeric(loc[s1, ])
    t = as.numeric(loc[s2, ])
    out = Kst(s, t, smothness, b1, b2)
    matern.mat[s1, s2] = out
  }
  matern.mat[lower.tri(matern.mat)] = t(matern.mat)[lower.tri(matern.mat)]
  #   if(!is.positive.definite(matern.mat))
  #     matern.mat = matrix(nearPD(matern.mat, corr = TRUE, keepDiag = TRUE)$mat, nrow(matern.mat), ncol(matern.mat))
  matern.mat
}
