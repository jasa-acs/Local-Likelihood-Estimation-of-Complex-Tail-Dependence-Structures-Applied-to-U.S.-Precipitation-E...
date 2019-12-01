### =================================================================== ###
### Auxiliary functions for the U.S. precipitation extremes application ###
### Daniela Castro-Camilo                                               ###
### Email: daniela.castro.camilo@gmail.com                              ###
### =================================================================== ###

#################
### Functions ###
#################
# F1: computed the univarite cdf for the stationary factor copula model.
# dummy: auxiliar function to compute the inverse of F1.
# F1inv: the inverse of F1.
# mi_condMVN : Computes conditional mean and covariance to be used in logCopula.R
# u.np: Function to non-parametrically transform to uniform margins

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

u.np = function(x) {
  rank(x, ties.method = "random", na.last = "keep")/(length(x[!is.na(x)]) + 1)
}

