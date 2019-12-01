### ==================================================================================================== ###
### Plot of empirical quantiles of five day-cumulative winter precipitation data (Figure 6 in the paper) ###                              
### Daniela Castro-Camilo                                                                                ###
### Email: daniela.castro.camilo@gmail.com                                                               ###
### ==================================================================================================== ###
print('Example usage: plot 90% quantile in log scale')
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

##########################
### Data and locations ###
##########################
load('DataApplication/Data/datamat5days1218.Rdata')
coord = read.table('DataApplication/Data/locations1218.txt', header = T)
datamat.90 = apply(datamat, 2, quantile, probs = 0.90, na.rm = TRUE)

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
nsbrks = c(30, 40, 50)
ewlbls <- unlist(lapply(ewbrks, function(x) paste(abs(x), "°W")))
nslbls <- unlist(lapply(nsbrks, function(x) paste(x, "°N")))
mycols = rev(tim.colors(25))

######################################
### Plot 90% quantile in log scale ###
######################################
dat = data.frame(long = c(10, coord$long,10), lat =c(10, coord$lat,10), value = c(1.8, log(datamat.90),6.4)); range(dat$value)

gg <- ggmap(map, extent = "device", legend = "bottomright", padding = 0.03) 
gg <- gg + geom_point(aes(x = dat$long, y = dat$lat, colour = dat$value), data = dat, alpha = 1, na.rm = T, shape = 16, size = 4)
gg <- gg + scale_color_gradientn(colours = mycols)
gg <- gg + geom_map(data = us, map = us, aes(x = long, y = lat, map_id = region), color = 1, alpha = 0.1, size = .7)
gg <- gg + labs(colour = "90% quantile")
gg <- gg + scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0), limits = c(-126, -66))
gg <- gg + scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0), limits = c(24, 51))
gg <- gg + theme(axis.text = element_text(size= 18))
gg <- gg + theme(axis.ticks = element_line(size = 1))
gg <- gg + theme(axis.ticks.length = unit(.1, "cm"))
gg <- gg + theme(legend.text = element_text(size = 15))
gg <- gg + theme(legend.key.width = unit(0.7, "cm"))
gg <- gg + theme(legend.position=c(1, 0.01))
gg <- gg + theme(legend.title = element_text(face = "bold", size = 20))
ggsave(gg, filename = 'Results/emp90q.pdf', width = 12, height = 7)

