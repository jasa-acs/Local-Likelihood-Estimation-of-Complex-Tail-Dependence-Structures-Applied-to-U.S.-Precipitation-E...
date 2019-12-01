### ========================================================================= ###
### Code to put 5-days cumulative winter precipitation data in a N x n matrix ###
### Daniela Castro-Camilo                                                     ###
### Email: daniela.castro.camilo@gmail.com                                    ###
### ========================================================================= ###
warning('Need to unzip all files in PreProcessing/RawData and then run DataApplication/PreProcessing/1_Preprocessing.R. This example usage will only produce results for state01_AL and it will finish with the error <cannot open compressed file DataApplication/PreProcessing/States/state02_AZ_6.Rdata, probable reason No such file or directory>.')

#################################################
###                                           ###
### Merge all stations and select winter data ###
###                                           ###
#################################################
require(rje)

states = c('01_AL', '02_AZ', '03_AR', '04_CA','05_CO', '06_CT', '07_DE', '08_FL', '09_GA', '10_ID', 
           '11_IL','12_IN', '13_IA', '14_KS', '15_KY', '16_LA', '17_ME', '18_MD', '19_MA','20_MI', 
           '21_MN', '22_MS', '23_MO', '24_MT', '25_NE', '26_NV', '27_NH', '28_NJ', '29_NM', '30_NY',
           '31_NC', '32_ND', '33_OH', '34_OK', '35_OR', '36_PA', '37_RI', '38_SC', '39_SD', '40_TN',
           '41_TX', '42_UT', '43_VT', '44_VA', '45_WA', '46_WV', '47_WI', '48_WY')

datamat = NULL
agnos = 1900:2014
ddays = 5

for(iter in 1:length(states)){
  printPercentage(iter, length(states))
  i = states[iter]
  temp = paste0('DataApplication/PreProcessing/States/state',i,'_6','.Rdata')
  load(temp)
  data1 = dat1[dat1[ , 2] >= agnos[1], ]
  data2 = data1[data1[ , 3] == 1 | data1[ , 3] == 2 | data1[ , 3] == 3 | data1[ , 3] == 12, ] # december to march
  station = unique(data2[, 1])
  data02 = NULL
  for(j in 1:length(station)){# for each station
    data3 = data2[ data2[ , 1] == station[j], ]
    years = unique(data3[, 2])
    if(length(years) < length(agnos)){
      aux1 = agnos[-match(years, agnos)] # missing years
      data3. = data3
      for(m in 1:length(aux1)){
        aux2 = rbind(c(station[j], aux1[m], 1, rep(-9999, 31)), c(station[j], aux1[m], 2, rep(-9999, 31)),
                     c(station[j], aux1[m], 3, rep(-9999, 31)), c(station[j], aux1[m], 12, rep(-9999, 31)))
        colnames(aux2) = colnames(data3.)
        data3. = rbind(data3., aux2)
      }
      data3. = data3.[order(data3.[, 2]), ]
      data3 = data3.  
    }
    data01 = fechas = NULL
    years = unique(data3[, 2])
    for(k in 1:length(years)){# for each year, take the 5 days cumulative rain
      data4 = data3[data3[ , 2] == years[k], ]
      if(nrow(data4) < 4){ # if there are missing months
        meses = c(1,2,3,12)
        tmp1 = meses[-match(data4[ , 3], meses)] # missing months
        data4. = data4
        for(l in 1:length(tmp1))
          data4. = rbind(data4., c(station[j] ,years[k], tmp1[l], rep(-9999, 31)))
        data4. = data4.[order(data4.[ , 3]), ]
        data4 = data4.
      }
      # Cleaning: remove feb 29 from leap years, select only winter days from december and march
      data4[data4[ , 3] == 2, 32] = NA; data4[data4[ , 3] == 2, 33] = NA; data4[data4[ , 3] == 2, 34] = NA
      data4[data4[ , 3] == 3, 24:34] = NA; data4[data4[ , 3] == 12, 4:23] = NA
      # 5-days cumulative rain over winter, for a given year in a given station
      data5 = data4[, -(1:3)]
      xmes = data4[ , 3]
      Xmes = NULL
      for(i in 1:length(xmes))
        Xmes = c(Xmes, rep(xmes[i], ncol(data5)))
      x = as.numeric(unlist(t(data5)))
      Xmes = Xmes[!is.na(x)]; x = x[!is.na(x)]
      Xmes[x == -9999] = NA ; x[x == -9999] = NA 
      x5 = as.numeric(unname(tapply(x, (seq_along(x)-1) %/% ddays, sum))) #
      fechas = c(fechas, rep(years[k], length(x5)))
      data01 = c(data01, x5)
    }
    data02 = cbind(data02, data01)
  } 
  datamat = cbind(datamat, data02)
  save(datamat, file = "Results/datamat5days1218.Rdata")
}
print('Done. The output file is a 2070x1218 matrix containing the required data to fit the model.')


