### ============================================== ###
### Code to remove observations with quality flags ###                              
### Daniela Castro-Camilo                          ###
### Email: daniela.castro.camilo@gmail.com         ###
### ============================================== ###
warning('Need to unzip all files in DataApplication/PreProcessing/RawData. This example usage will only produce results for state01_AL and it will finish with the error <cannot open file DataApplication/PreProcessing/RawData/state02_AZ.txt: No such file or directory>')

# -------------------------------------------------------------------------------------------------
# UNITED STATES HISTORICAL CLIMATOLOGY NETWORK (USHCN) Daily Dataset
# M.J. Menne, C.N. Williams, Jr., and R.S. Vose
# National Climatic Data Center, National Oceanic and Atmospheric Administration

# These files comprise CDIAC's most current version of USHCN daily data.
# Data are available as 48 individual ASCII state files 
# and as five netcdf files, both with data through 2014.
# The netcdf files are organized by specific data type (i.e., ushcn_datatype.nc.gz;
# see directory contents below).
# More info: http://cdiac.ornl.gov/epubs/ndp/ushcn/ushcn.html
#------------------------------------------------------------------------------------------

ns = c('01_AL', '02_AZ', '03_AR', '04_CA','05_CO', '06_CT', '07_DE', '08_FL', '09_GA', '10_ID', 
       '11_IL','12_IN', '13_IA', '14_KS', '15_KY', '16_LA', '17_ME', '18_MD', '19_MA','20_MI', 
       '21_MN', '22_MS', '23_MO', '24_MT', '25_NE', '26_NV', '27_NH', '28_NJ', '29_NM', '30_NY',
       '31_NC', '32_ND', '33_OH', '34_OK', '35_OR', '36_PA', '37_RI', '38_SC', '39_SD', '40_TN',
       '41_TX', '42_UT', '43_VT', '44_VA', '45_WA', '46_WV', '47_WI', '48_WY')
names = NULL
for(i in 1:length(ns))
  names[i] = paste('state', ns[i], '.txt',sep = '')

for(i in 1:length(names)){
  print(names[i])
  wd. = paste0('DataApplication/PreProcessing/RawData/', names[i])
  dat = read.fwf(wd.,widths = c(6,4,2,4,rep(c(5,1,1,1),31)),na.strings = " ")
  name1 = paste('state', ns[i], sep = '')
  # Only precipitation data
  dat1 = dat[dat[ ,4]=="PRCP", ]
  # the first step was identifying all observations in the 1218 USHCN station
  # records having any of the quality flag assigments described on the 
  # quality control page of the GHCN documentation (see also Durre et al. 2010). 
  # These flags indicate that the accuracy and reliability of data are 
  # questionable, so these observations were set to the missing indicator 
  # ("-999") matching the one already used in the database for actual missing 
  # observations.
  # V7, V11, V15, V19,... = quality flag
  i = seq(0,((127-7)/4),1) # all quality flags
  temp1 = dat1[ ,(7 + 4 * i)]
  # temp1[1,1] = temp1[3,3] = temp1[4,5] = temp1[1,10] = 1 # to check if works
  temp2 = dat1[ ,(5 + 4 * i)]
  temp2[!is.na(temp1)] = -9999
  dat1[ ,(5 + 4 * i)] = temp2
  
  # Remove all the columns that we don't need:
  dat1 = dat1[ ,c(1,2,3,(5 + 4 * i))]
  # wd1 = paste0('PreProcessing/States/',name1,'_6', '.Rdata')
  wd1 = paste0('Results/',name1,'_6', '.Rdata')
  save(dat1, file = wd1)
}
print('Done.')



