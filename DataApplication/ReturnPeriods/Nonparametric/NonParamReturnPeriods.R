### ============================================================================================================ ###
### Compute non-parametric return periods associated with catastrophic events over selected stations in 5 states ###
### Daniela Castro-Camilo                                                                                        ###
### Email: daniela.castro.camilo@gmail.com                                                                       ###
### ============================================================================================================ ###

#################################################################
### Libraries, auxiliary codes, and data at selected stations ###
#################################################################
load('DataApplication/Data/datamat5days1218.Rdata')
row.na = function(r) any(is.na(r))
df = read.table('DataApplication/ReturnPeriods/Parametric/param_selectedstations_nu=0.5.txt', header = T)
datamat = datamat[, df$ID] # data over the stations

##################################################################
### Function to compute non-parametric return period estimates ###
##################################################################
get.np.rp = function(data, umax){
  is.na = apply(data, 1, row.na)
  xmat1 = data[!is.na, ]
  u = as.numeric(apply(xmat1, 2, quantile, probs = umax))
  
  N = nrow(xmat1)
  count = rep(NA, N)
  for(i in 1:N)
    count[i] = as.numeric(xmat1[i,1]) > u[1] & as.numeric(xmat1[i,2]) > u[2] & as.numeric(xmat1[i,3]) > u[3]
  
  # Point estimate
  phat = mean(count)
  if(phat == 0)
    return('phat = 0')
  else{
    rphat = ((1/phat) * 5)/365 # return period in years
    return(list(rphat = rphat))
  }
}

############################################
### Function to compute sd via bootstrap ###
############################################
get.np.rp.boot = function(b, data, umax){
  is.na = apply(data, 1, row.na)
  data = data[!is.na, ]
  id.boot = sample(1:nrow(data), nrow(data), replace = T)
  data.boot = data[id.boot, ]
  out = get.np.rp(data.boot, umax)
  out$rphat
}

##############################################
### Non-parametric return period estimates ###
##############################################
print('Computing non-parametric return period estimates and bootstrap-based confidence intervals for u=0.94')
umax = 0.94
B = 1000

np.rp = function(s, data, umax, B){
  data = data[, s]
  out = get.np.rp(data, umax)
  rpboot = sapply(1:B, get.np.rp.boot, data = data, umax = umax)
  out$ll.boot = out$rphat - 2 * sd(rpboot)
  out$ul.boot = out$rphat + 2 * sd(rpboot)
  out
}

# Loussiana
out.LA = np.rp(4:6, datamat, umax, B)
# Mississippi
out.MS = np.rp(7:9, datamat, umax, B)
# Kentucky
out.KY = np.rp(10:12, datamat, umax, B)
# Florida
out.FL = np.rp(16:18, datamat, umax, B)
# Tennessee
out.TN = np.rp(19:21, datamat, umax, B)


rp = data.frame(matrix(unlist(c(out.LA, out.MS, out.KY, out.FL, out.TN)), ncol = 3, byrow = T))
colnames(rp) = c('rp.hat', 'lb.boot', 'up.boot')
rp = round(rp,2)
rp$State = c('LA', 'MS', 'KY', 'FL', 'TN')
print('Non-parametric return periods')
rp


