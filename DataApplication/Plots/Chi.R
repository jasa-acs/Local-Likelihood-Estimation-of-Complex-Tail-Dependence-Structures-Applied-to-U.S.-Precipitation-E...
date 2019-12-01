### ============================================================== ###
### Compute and plot estimates of chi_h(u) (Figure 7 in the paper) ###                              
### Daniela Castro-Camilo                                          ###
### Email: daniela.castro.camilo@gmail.com                         ###
### ============================================================== ###
print('Example usage: compute and plot estimates of chi_h(u) for h = 40, u = 0.95 (2nd row, 1st column in Figure 7 in the paper)')
########################################
### Ignore if not using docker image ###
########################################
args <- commandArgs(TRUE)
for (arg in args) {
  eval(parse(text = arg))
}
rm(arg, args)

#####################################
### Libraries and auxiliary codes ###
#####################################
library(ggmap)
library(fields)
library(ggplot2)
library(mvtnorm)
source('DataApplication/Fit/Tools.R')

#####################################################################################
### Function to compute chi_h(u) for given rate, range, and smoothness parameters ###
#####################################################################################
chi.m_h = function(h, u, lbda, rho, nu){
  dist.mat = matrix(c(0,h,h,0), nrow = 2, ncol = 2)
  M = Matern(dist.mat, range = rho, smoothness = nu)
  r = M[1,2]
  F2 = function(z, lbda, r){
    p1 = pmvnorm(lower = -Inf, upper = c(z,z), mean = rep(0,2), sigma = matrix(c(1, r, r, 1), ncol =2))[1] 
    S0 = matrix(c(2 * (1 - r), -(1 - r), -(1 - r), 1), nrow = 2)
    z0 = c(z - r*z, 0)
    mu0 = c((z - lbda) * (1 - r), (lbda - z))
    p2 = pmvnorm(lower = -Inf, upper = z0, mean = mu0, sigma = S0)[1]
    return(p1 - 2 * exp(lbda^2/2 - lbda * z) * p2)
  }
  z = F1inv(u, lbda = lbda)
  cu = F2(z, lbda = lbda, r = r)
  chiu = (1 - 2*u + cu)/(1 - u)
  chiu
}

###############################################
### Computing chi_h(u) for h = 40, u = 0.95 ###
###############################################
param = read.table('DataApplication/Plots/param_nu=0.5.txt', header = T)
nu = 0.5
h = 40
u = 0.95

chi = numeric(nrow(param))
for(i in 1:nrow(param)){
  rate = param$param1[i]
  range = param$param2[i]
  chi[i] = chi.m_h(h, u, rate, range, nu = nu)
}

########################
### Set map features ###
########################
warning('Introduce your api key to enable Google services in R')
key = api_key # replace with your api key to enable Google services in R. Ignore if using docker image.
register_google(key = key)
map <- get_map(location = c(-126, 27, -60, 48), zoom = 3, maptype = "satellite", source = 'google', color = 'color')
us <- map_data("state")

###########################################
### create the breaks and label vectors ###
###########################################
ewbrks = c(-120, -110, -100, -90, -80)
nsbrks = c(20, 40, 50)
ewlbls <- unlist(lapply(ewbrks, function(x) paste(abs(x), "°W")))
nslbls <- unlist(lapply(nsbrks, function(x) paste(x, "°N")))
mycols = rev(tim.colors(25))

##########################################
### Plot chi_h(u) for h = 40, u = 0.95 ###
##########################################
dat = data.frame('lon' = param$lon, 'lat' = param$lat, 'value' = chi)
# To be able to have a scale from 0 to 1 included:
dat = rbind(dat, c(10, 10, 0))
dat = rbind(dat, c(10, 10, 1))

gg <- ggmap(map, extent = "device", legend = "bottomright", padding = 0.03) 
gg <- gg + geom_point(aes(x = dat$lon, y = dat$lat, colour = dat$value), data = dat, alpha = 1, na.rm = T, shape = 15, size = 8)
gg <- gg + scale_color_gradientn(colours = tim.colors(25))
gg <- gg + geom_map(data = us, map = us, aes(x = long, y = lat, map_id = region), color = 1, alpha = 0.1, size = .7)
gg <- gg + labs(colour = eval(parse(text = paste0('expression(hat(chi)[', 40, '](', 0.95, '))'))))
gg <- gg + scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0), limits = c(-126, -66))
gg <- gg + scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0), limits = c(24, 51))
gg <- gg + theme(axis.text = element_text(size= 18))
gg <- gg + theme(axis.ticks = element_line(size = 1))
gg <- gg + theme(axis.ticks.length = unit(.1, "cm"))
gg <- gg + theme(legend.text = element_text(size = 15))
gg <- gg + theme(legend.key.width = unit(0.7, "cm"))
gg <- gg + theme(legend.position = c(1, 0.01))
gg <- gg + theme(legend.title = element_text(face = "bold", size = 20))
ggsave(gg, filename = 'Results/Chi_h40_u095.pdf', width = 12, height = 7)


