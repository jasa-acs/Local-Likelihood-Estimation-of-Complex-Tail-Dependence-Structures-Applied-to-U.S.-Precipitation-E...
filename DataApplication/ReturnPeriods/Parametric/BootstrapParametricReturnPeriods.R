### ============================================================================================================================================== ###
### Bootstrap to compute sd of parametric (simulation-based) return periods associated with catastrophic events over selected stations in 5 states ###
### Includes simulation of N replicates from our model and transformation to uniform margins                                                       ###
### Daniela Castro-Camilo                                                                                                                          ###
### Email: daniela.castro.camilo@gmail.com                                                                                                         ###
### ============================================================================================================================================== ###

#################################################################
### Libraries, auxiliary codes, and data at selected stations ###
#################################################################
library(rje)
library(mvtnorm)
library(matrixcalc)
library(fields)
library(ismev)
library(evd)
library(parallel)
load("DataApplication/Data/datamat5days1218.Rdata")
df = read.table('DataApplication/ReturnPeriods/Parametric/param_selectedstations_nu=0.5.txt', header = T)
coord = read.table('DataApplication/Data/locations1218.txt', header = T) ; coord = coord[ , 2:3]
nu = 0.5

#########################################################################################
### Auxiliary: Simulating a Gaussian Process with non stationary correlation function ###
#########################################################################################
Mnu = function(arg, nu)
  abs(arg) ^ nu * besselK(arg, nu)

rM.st = function(locs, rho1, rho2, nu){
  s = locs[1, ]; t = locs[2, ]
  if((s[1] == t[1]) & (s[2] == t[2]))
    return(1)
  tmp = rep(NA, 3)
  d = rdist.earth(locs, miles = F)[1,2]
  if (any(d < 0)) 
    stop("distance argument must be nonnegative")
  d[d == 0] <- 1e-10
  d = d/(2*sqrt(nu))
  con = (2^(nu - 2)) * gamma(nu)
  tmp[1] = 1/con
  den = rho1^2 + rho2^2
  tmp[2] = (rho1*rho2)/den
  arg = (2 * sqrt(2 * nu) * d)/sqrt(den)
  tmp[3] = Mnu(arg, nu)
  return(prod(tmp))
}

################################################################################################################
### Simulating N replications from our model and transforming to uniform margins, for every bootstrap sample ###
################################################################################################################
print('Example usage: simulating N = 5e04 replications from our model and transforming to uniform margins, for each of B=10 bootstrap sample. In the paper we use N = 5e05 and B=300 samples.')
D = nrow(df)
lbda = df$param1
M = matrix(NA, D, D); diag(M) = (lbda)^2
for(i in 1:D){
  for(j in 1:D){
    if(i < j) M[i,j] = rM.st(coord[c(i,j), ], rho1 = df$param2[i], rho2 = df$param2[j], nu)
  }
}

M[lower.tri(M)] = t(M)[lower.tri(M)]
# is.symmetric.matrix(M)
if(!is.positive.definite(M)){
  M. = nearPD(M, corr = TRUE)
  M = as.matrix(M.$mat)
  chol.M = chol(M)
}else{
  chol.M = chol(M)
}

N = 5e04 # N=5e05 in the paper
B = 10 # B=300 in the paper

bootsim = function(b, N, M, lbda){
  set.seed(b) 
  Z = rmvnorm(N, sigma = M)
  set.seed(2*b)
  V = rexp(N, rate = lbda)
  sim = Z + V
}

F1 = function(w, lbda){
  pnorm(w) - pnorm(w - lbda) * exp(lbda^2/2 - lbda * w)
}

bootsimu = function(b, N, D, lbda){
  tmp = sim[[b]]
  tmp.u = matrix(NA, N, D)
  for(i in 1:N){
    w = tmp[i, ]
    tmp.u[i, ] = mapply(F1, w, lbda)
  }
  tmp.u
}

sim = Data.u = list()
mccores = detectCores() 
njobs = ceiling(B/mccores)
for(i in 1:njobs){
  resto = B%%mccores
  tmp = if ( i == njobs && resto > 0 ) resto else mccores
  x = (mccores)*(i-1) + c(1:(tmp))
  mccores = tmp
  
  sim = c(sim, mclapply(x, bootsim, N = N, M = M, lbda = lbda, mc.cores = mccores))
}

for(i in 1:(B/mccores)){
  resto = B%%mccores
  tmp = if ( i == njobs && resto > 0 ) resto else mccores
  x = (mccores)*(i-1) + c(1:(tmp))
  mccores = tmp
  
  Data.u = c(Data.u, mclapply(x, bootsimu, N = N, D = D, lbda = lbda, mc.cores = mccores))
}

#####################################################################
### Parametric return period estimates for every bootstrap sample ###
#####################################################################
print('Computing parametric return period estimates for every bootstrap sample using u=0.94')
p.rp = function(s, data, data.u, u.max){
  # s: station
  # data: simualted observations
  # data.u: simulated observations in uniform scale
  # u.max: threshold, e.g., 0.995
  
  # Auxiliar function: GPD fit to stabilize resuts
  qtilde = function(u, xdat, prob = 0.95){
    # Obtain the u%-quantile of the GPD
    xdat = xdat[!is.na(xdat)]
    th = quantile(xdat, prob)
    fit = gpd.fit(xdat = xdat, threshold = th, show = FALSE)
    s = fit$mle[1]; xi = fit$mle[2]
    qgpd(p = u, scale = s, shape = xi)
  }
  
  # Obtain the thresholds in the scale of the data
  w = NULL
  for(i in s){
    xdat = data[ ,i]; xdat = xdat[!is.na(xdat)]
    out = qtilde(u.max, xdat)
    w = c(w, out) # thresholds
  }
  
  U = data.u[, s]
  rowmin = apply(U, 1, min)
  pr = mean(rowmin > u.max)
  rp = ((1/pr) * 5)/365 # return period in years
  list(pr = pr, rp = rp, thresholds = w)
}

p.rp.boot = function(b, s, data, Data.u, u.max){
  data.u = Data.u[[b]]
  out = p.rp(s, data, data.u, u.max)
  out$rp
}
  
umax = 0.94
# Loussiana
out.LA = sapply(1:B, p.rp.boot, s = 1:3, data = datamat, Data.u = Data.u, u.max = umax)
# Mississippi
out.MS = sapply(1:B, p.rp.boot, s = 7:9, data = datamat, Data.u = Data.u, u.max = umax)
# Kentucky
out.KY = sapply(1:B, p.rp.boot, s = 10:12, data = datamat, Data.u = Data.u, u.max = umax)
# Florida
out.FL = sapply(1:B, p.rp.boot, s = 16:18, data = datamat, Data.u = Data.u, u.max = umax)
# Tennessee
out.TN = sapply(1:B, p.rp.boot, s = 19:21, data = datamat, Data.u = Data.u, u.max = umax)

###########################################
### Bootstrap-based standard deviations ###
###########################################
print('Bootstrap-based standard deviations for parametric return period estimates using u=0.94')
sd.rp = c(sd(out.LA), sd(out.MS), sd(out.KY), sd(out.FL), sd(out.TN))
sd.rp = data.frame(sd.rp = sd.rp, State = c('LA', 'MS', 'KY', 'FL', 'TN'))
sd.rp

