### ========================================================================================================== ###
### Boxplots of sd for estimates chi_h(u) for  h = 20, 40, 80, 160 and u = 0.95, 0.98 (Figures 8 in the paper) ###                              
### Daniela Castro-Camilo                                                                                      ###
### Email: daniela.castro.camilo@gmail.com                                                                     ###
### ========================================================================================================== ###
print('Example usage: boxplots of sd for estimates chi_h(u) for  h = 20, 40, 80, 160 and u = 0.95, 0.98 (Figures 8 in the paper)')
########################################
### Ignore if not using docker image ###
########################################
args <- commandArgs(TRUE)
for (arg in args) {
  eval(parse(text = arg))
}
rm(arg, args)

####################################################
### Libraries, auxiliary codes, and data to plot ###
####################################################
library(ggmap)
library(fields)
library(ggplot2)
library(mvtnorm)
param = read.table('DataApplication/Plots/sdchi_nu=0.5.txt', header = T)

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
mycols = "royalblue4"

#############################################################
### Plot chi_h(u) for h in {20, 40, 80, 160} and u = 0.95 ###
#############################################################
n = nrow(param)
dat = data.frame('label' = c(rep("1", n), rep("2", n), rep("3", n), rep("4", n)), 'value' = c(param$sd_chi_u95_h20, param$sd_chi_u95_h40, param$sd_chi_u95_h80, param$sd_chi_u95_h160))
gg <- ggplot(data = dat, aes(x = label, y = value))
gg <- gg + geom_boxplot(fill = mycols, colour = "#1F3552", alpha = 0.7, outlier.colour = "#1F3552", outlier.shape = 20)
gg <- gg + scale_x_discrete(name = "Distance (km)", labels = c("h = 20", "h = 40", "h = 80", "h = 160")) + scale_y_continuous(name = "")
gg <- gg + ylim(c(0, 0.18)) + ylab('')
gg <- gg + theme(axis.text = element_text(size = 20), 
                 axis.title = element_text(size = 20), 
                 plot.title = element_text(hjust = 0.5, size = 30))
gg <- gg + ggtitle('u=0.95')
ggsave(gg, filename = 'Results/sdChi_u095.pdf', width = 9, height = 9)


#############################################################
### Plot chi_h(u) for h in {20, 40, 80, 160} and u = 0.98 ###
#############################################################
dat = data.frame('label' = c(rep("1", n), rep("2", n), rep("3", n), rep("4", n)), 'value' = c(param$sd_chi_u98_h20, param$sd_chi_u98_h40, param$sd_chi_u98_h80, param$sd_chi_u98_h160))
gg <- ggplot(data = dat, aes(x = label, y = value))
gg <- gg + geom_boxplot(fill = mycols, colour = "#1F3552", alpha = 0.7, outlier.colour = "#1F3552", outlier.shape = 20)
gg <- gg + scale_x_discrete(name = "Distance (km)", labels = c("h = 20", "h = 40", "h = 80", "h = 160")) + scale_y_continuous(name = "")
gg <- gg + ylim(c(0, 0.18)) + ylab('')
gg <- gg + theme(axis.text = element_text(size = 20), 
                 axis.title = element_text(size = 20), 
                 plot.title = element_text(hjust = 0.5, size = 30))
gg <- gg + ggtitle('u=0.98')
ggsave(gg, filename = 'Results/sdChi_u098.pdf', width = 9, height = 9)

