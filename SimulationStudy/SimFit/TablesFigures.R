### ============================================================================= ###
### Guidelines to reproduce Tables and Figure in the simulation study (Section 4) ###
### Daniela Castro-Camilo                                                         ###
### Email: daniela.castro.camilo@gmail.com                                        ###
### ============================================================================= ###

#####################################
### Libraries and auxiliary codes ###
#####################################
library(plot3D)
library(fda)
source('SimulationStudy/SimFit/Tools.R')
# Function to compute bivariate chi measure for true rate and range parameters
chi = function(s, t, b1, b2, nu){
  lbda.s = lambda.st(s, nu, b1, b2)
  lbda.t = lambda.st(t, nu, b1, b2)
  r = Kst(s, t, nu, b1, b2)
  gamma = lbda.s ^ 2 + lbda.t ^ 2 - 2 * lbda.s * lbda.t * r
  2 - 2 * pnorm(gamma^(1/2)/2)
}
# Grid to plot results
grid.length = 100
xgrid = ygrid = seq(1, 10, length.out = floor(sqrt(grid.length)))
grid = expand.grid(xgrid, ygrid)

########################################################################################################################
### To reproduce Table 1, generate 1000 replicates and fit the model using fit.sim with the following configurations ###
########################################################################################################################
# Weakly non-stationary: b1 = 5; b2 = 4 with nu = 0.5, 1.5, and 2.5
# Midly non-stationary: b1 = 5; b2 = 1.5 with nu = 0.5, 1.5, and 2.5
# Strongly non-stationary: b1 = 5; b2 = .8 with nu = 0.5, 1.5, and 2.5

########################################################################################################################
### To reproduce Table 2, generate 1000 replicates and fit the model using fit.sim with the following configurations ###
########################################################################################################################
# In SimFit/ConsoredLocalLikelihood.R, replace line 35 by lbda.st[i] = lambda.st(as.numeric(coord[i, ]), smooth, b1 = 5, b2 = 0.8)
# Weakly non-stationary: b1 = 5; b2 = 4 with nu = 2.5 and neigh = 5, 10, 15, 20, 25
# Midly non-stationary: b1 = 5; b2 = 1.5 with nu = 2.5 and neigh = 5, 10, 15, 20, 25
# Strongly non-stationary: b1 = 5; b2 = .8 with nu = 2.5 and neigh = 5, 10, 15, 20, 25


###############################################################
### Figure 4: Plot of true surfaces (mildly non-stationary) ###
###############################################################
print('Plot of true surfaces for the mildly non-stationary')
# True values of rate, range, and bivariate measure chi
b1 = 5; b2 = 1.5; nu = 2.5
x = y = seq(1, 10, length.out = 100)
sx = grid[45,1]; sy = grid[45,2]; s = c(sx, sy)
lbda.s = lambda.st(s, nu, b1, b2)
rho.s = rho(s, b1, b2)
delta.true  = lambda.true = chi.true =  matrix(NA, length(x), length(y))
for(i in 1:length(x)){
  for(j in 1:length(y)){
    delta.true[i, j] = rho(c(x[i], y[j]), b1, b2)
    lambda.true[i, j] = lambda.st(c(x[i], y[j]), nu, b1, b2)
    chi.true[i,j] = chi(s, c(x[i], y[j]), b1, b2, nu)
  }
}

# Plot rate lambda
pdf('Results/TrueLogRateMildly.pdf', width = 10,height = 7)
persp3D(x = x, y = y, z = log(lambda.true), scale = T, expand = 1, bty = "g", phi = 0,
        zlim = c(-1.5,1.5), col = "royalblue4", lighting = T, ltheta = 50, along = "y", space = 0.7, ticktype = "detailed",
        d = 2, curtain = F, ylab = 'Sy', xlab = 'Sx', zlab = '', colkey = F, add = F, theta = 50, axes = T)
mtext(expression(log(lambda[s])), side = 2, cex = 1.7, line = -7)
dev.off()

# Plot range delta
pdf('Results/TrueRangeMildly.pdf', width = 10,height = 7)
persp3D(x = x, y = y, z = delta.true, scale = T, expand = 1, bty = "g", phi = 0,
        zlim = c(0.5, 1.6), col = "royalblue4", lighting = T, ltheta = 50, along = "y", space = 0.7, ticktype = "detailed",
        d = 2, curtain = F, ylab = 'Sy', xlab = 'Sx', zlab = '', colkey = F, add = F, theta = 35, axes = T)
mtext(expression(delta[s]), side = 2, cex = 1.7, line = -7)
dev.off()

# Plot bivariate measure chi (s = (5,5))
pdf('Results/TrueChiMildly_s=(5,5).pdf', width = 10,height = 7)
persp3D(x = x, y = y, z = chi.true, scale = T, expand = 1, bty = "g", phi = 0,
        zlim = c(0, 1), col = "royalblue4", lighting = T, ltheta = 50, along = "y", space = 0.7, ticktype = "detailed",
        d = 2, curtain = F, ylab = 'Sy', xlab = 'Sx', zlab = '', colkey = F, add = F, theta = 50, axes = T)
mtext(expression(chi[12]), side = 2, cex = 1.7, line = -7)
dev.off()

#################################################
### Figure 5: functional and surface boxplots ###
#################################################
print('Ploting functional and surface boxplots')
# Functional boxplot using a sample of average estimated range at the 10 locations
load('SimulationStudy/SimFit/avg_range.Rdata')
pdf('Results/fbplotDelta.pdf', width = 10,height = 7)
fbp = fbplot(avr.range, xlim = c(1, 10), color = c('mediumpurple'), ylab = expression(rho[s[x]]), method = 'Both', 
             plot = T, factor = 1.3)
dev.off()

# Each of the components of the surface boxplot for the bivariate chi measure
load('SimulationStudy/SimFit/chi_centerL.Rdata')
load('SimulationStudy/SimFit/chi_centerU.Rdata')
load('SimulationStudy/SimFit/chi_median.Rdata')

pdf('Results/sbplotChi.pdf', width = 10,height = 7)
par(pty = "m", cex.lab = 2, cex.axis = 1.5, mar = c(4, 3, 2, 0), oma = c(1, 1, 1, 1))
persp3D(x = xgrid, y = ygrid, z = chi.centerL, scale = T, expand = 1, bty = "g", phi = 0,
        zlim = c(0,1), col = 'darkturquoise', lighting = T, ltheta = 50, along = "y", space = 0.7, 
        ticktype = "detailed", d = 2, curtain = F, ylab = 'Sy', xlab = 'Sx', zlab = '', colkey = F, theta = 50)

persp3D(x = xgrid, y = ygrid, z = chi.median, scale = T, expand = 1, bty = "g", phi = 0,
        col = 'chartreuse2', lighting = T, ltheta = 50, along = "y", space = 0.7, ticktype = "detailed",
        d = 2, curtain = F, colkey = F, add = T, theta = 50)

persp3D(x = xgrid, y = ygrid, z = chi.centerU, scale = T, expand = 1, bty = "g", phi = 0,
        col = 'royalblue4', lighting = T, ltheta = 50, along = "y", space = 0.7, ticktype = "detailed",
        d = 2, curtain = F, colkey = F, add = T, theta = 50)

mtext(expression(hat(chi)[12]), side = 2, cex = 1.5, line = 2)
mtext(expression(s[2] ~" = (2,3)"^T), side = 3, cex = 1.7, line = 0)
dev.off()
