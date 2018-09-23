library(synergyfinder)
library(plyr)
library(openxlsx)
library(mgcv)
library(dplyr)
library(SpatialExtremes)
library(reshape2)
library(drc)
setwd("C:/Users/Localadmin_jtang/Dropbox/drugcomb/code")

source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/reshape2.R')
options(scipen = 999) # suppress the scientific notation

# READ RESPONSE DATA
rm(list=ls(all=TRUE))
data = read.xlsx("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/database/response.xlsx")
own_rank = function(x){
  x_unique <- unique(x)
  x_ranks <- rank(x_unique)
  x_ranks[match(x,x_unique)]}

# add row and column numbers
data2 = ddply(data,c("cell_line_name","drug_row","drug_col","block_id"),transform, Row=own_rank(conc_r), Col=own_rank(conc_c))
names(data2)[2] = "BlockID"
names(data2)[3] = "ConcRow"
names(data2)[4] = "ConcCol"
names(data2)[5] = "Response"
names(data2)[6] = "DrugRow"
names(data2)[7] = "DrugCol"
names(data2)[8] = "ConcRowUnit"
names(data2)[9] = "ConcColUnit"

# response
data_response <- ReshapeData2(data2, data.type = "inhibition")

combo_response = list()
single_response = list()
for(i in unique(data$block_id)){
  cat('\r',i)
  flush.console()
  index = which(data_response$drug.pairs$blockIDs == i)
  
  # surface krigging
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_response$dose.response.mats[[index]], 1)
    }, 
    error = function(cond){
      print(i);
      mat1 = data_response$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1, 1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'response'))
  mat2$block_id = data_response$drug.pairs$blockIDs[index]
  mat2$drug_row = data_response$drug.pairs$drug.row[index]
  mat2$drug_col = data_response$drug.pairs$drug.col[index]
  mat2$conc_r_unit = data_response$drug.pairs$concRUnit[index]
  mat2$conc_c_unit = data_response$drug.pairs$concCUnit[index]
  combo_response[[i]] = mat2
  
  # single drug curve fitting
  mat3 = data_response$dose.response.mats[[index]]
  mat4 = setNames(melt(mat3), c('conc_r', 'conc_c', 'response'))
  
  mat5 = mat4[which(mat4$conc_r==0 | mat4$conc_c==0),]
  single.row = mat5[which(mat5$conc_c==0),]
  single.col = mat5[which(mat5$conc_r==0),]
  
  coef.row = tryCatch (
    { # in case error happens, return the score of NA
      tmp = drm(single.row$response ~ single.row$conc_r, fct=LL.4())
      tmp$coefficients
    }, 
    error = function(cond){
      print(i)
      return(c(NA,NA,NA,NA))
    }
  )
  
  coef.col = tryCatch (
    { # in case error happens, return the score of NA
      tmp = drm(single.col$response~single.col$conc_c, fct=LL.4())
      tmp$coefficients
    }, 
    error = function(cond){
      print(i)
      return(c(NA,NA,NA,NA))
    }
  )
  
  mat6 = data.frame(matrix(NA, nrow = 2, ncol = 9))
  mat6 = setNames(mat6, c("block_id","drug_row","drug_col","conc_r_unit","conc_c_unit","b","c","d","e"))

  mat6[1,] =  c(data_response$drug.pairs$blockIDs[index], data_response$drug.pairs$drug.row[index], NA,
                data_response$drug.pairs$concRUnit[index], NA, coef.row)                   
  mat6[2,] =  c(data_response$drug.pairs$blockIDs[index], NA, data_response$drug.pairs$drug.col[index], NA,
                data_response$drug.pairs$concCUnit[index], coef.col)  
  single_response[[i]] = mat6
    
  # fitted.row = coef.row[2] + (coef.row[3]-coef.row[2])/(1+exp(coef.row[1]*(log(x.row)-log(coef.row[4]))))
  # plot(curve.row)
  # points(x.row, fitted.row, col = "red")
  }
res_table_response = do.call(rbind, combo_response)
curve = do.call(rbind, single_response)
write.xlsx(curve,"curve.xlsx")

# zip score
colnames(data2)[5] = "NA"
colnames(data2)[10] = "Response"
data_zip<- ReshapeData2(data2, data.type = "inhibition")

res_zip = list()
for(i in unique(data$block_id)){
  index = which(data_zip$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_zip$dose.response.mats[[index]], 1)
    }, 
    error = function(cond){
      print(i);
      mat1 = data_zip$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1, 1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_zip'))
  mat2$block_id = data_zip$drug.pairs$blockIDs[i]
  res_zip[[i]] = mat2
}
res_table_zip = do.call(rbind, res_zip)

# bliss score
colnames(data2)[10] = "NA"
colnames(data2)[11] = "Response"
data_bliss<- ReshapeData2(data2, data.type = "inhibition")

res_bliss = list()
for(i in unique(data$block_id)){
  index = which(data_bliss$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_bliss$dose.response.mats[[index]], 1)
    }, 
    error = function(cond){
      print(i);
      mat1 = data_bliss$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1, 1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_bliss'))
  mat2$block_id = data_bliss$drug.pairs$blockIDs[i]
  res_bliss[[i]] = mat2
}
res_table_bliss = do.call(rbind, res_bliss)


# hsa score
colnames(data2)[11] = "NA"
colnames(data2)[13] = "Response"
data_hsa<- ReshapeData2(data2, data.type = "inhibition")

res_hsa = list()
for(i in unique(data$block_id)){
  index = which(data_hsa$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_hsa$dose.response.mats[[index]], 1)
    }, 
    error = function(cond){
      print(i);
      mat1 = data_hsa$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1, 1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_hsa'))
  mat2$block_id = data_hsa$drug.pairs$blockIDs[i]
  res_hsa[[i]] = mat2
}
res_table_hsa = do.call(rbind, res_hsa)


# loewe score
colnames(data2)[13] = "NA"
colnames(data2)[12] = "Response"
data_loewe<- ReshapeData2(data2, data.type = "inhibition")

res_loewe = list()
for(i in unique(data$block_id)){
  index = which(data_loewe$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_loewe$dose.response.mats[[index]], 1)
    }, 
    error = function(cond){
      print(i);
      mat1 = data_loewe$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1, 1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_loewe'))
  res_loewe[[i]] = mat2
}
res_table_loewe = do.call(rbind, res_loewe)

surface = cbind(res_table_response, res_table_zip, res_table_bliss, res_table_hsa, res_table_loewe)
write.csv(surface[,c(1:8,11,15,19,23)],"surface.csv")
save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/data_process.RData")


