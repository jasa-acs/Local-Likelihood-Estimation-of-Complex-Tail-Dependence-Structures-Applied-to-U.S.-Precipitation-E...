### ====================================== ###
### Generating block-Bootstrap samples     ###                              
### Daniela Castro-Camilo                  ###
### Email: daniela.castro.camilo@gmail.com ###
### ====================================== ###
print('Example usage: generating B=2 bootstrap samples. In the paper we use B=300.')
require(rje)
load('DataApplication/Data/datamat5days1218.Rdata')
B = 2 # B = 300 in the paper
N = nrow(datamat)
block.size = 6 # app one month in each block
n.blocks = floor(N/6)
id.blocks = matrix(NA, n.blocks, 3)
for(i in 1:n.blocks)
id.blocks[i, ] = c(i, (block.size*(i-1)+1), (block.size*i))
colnames(id.blocks) = c('block.number', "from", "to")

for(b in 1:B){
  printPercentage(b, B)
  new.datamat = NULL
  block.sample = sample(1:n.blocks, n.blocks, replace = T)
  row.sel = id.blocks[block.sample, -1]
  for(i in 1:n.blocks){
    temp = datamat[row.sel[i, 1]:row.sel[i, 2], ]
    new.datamat = rbind(new.datamat, temp)
  }
  name = paste0("datamatB", b, ".Rdata")
  # save(new.datamat, file = paste0("DataApplication/Fit/Bootstrap/BootstrapSamples/", name))
  save(new.datamat, file = paste0("Results/", name))
}

