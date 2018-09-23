library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
library(nleqslv)
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Oneil")
rm(list=ls(all=TRUE))

source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/FittingSingleDrug2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/CalculateSynergy2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ZIP2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Loewe2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Bliss2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/HSA2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/BaselineCorrectionSD2.R')
options(scipen = 999) # suppress the scientific notation

# READ RESPONSE DATA
response = read.xlsx("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Oneil/response_oneil_raw.xlsx")

# add row and column numbers
response = ddply(response, c("cell_line_name","drug_row","drug_col","block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))

scores = list()
m = length(unique(response$block_id))
options(show.error.messages = F)

for(i in 1:m){
  cat('\r',i)
  index = which(response$block_id == i)
  data.tmp = response[index,]
  data.tmp2 = ReshapeData2(data.tmp, data.type = "inhibition")
  
  # missing value imputation
  data.tmp3 = data.tmp2$dose.response.mats[[1]]
  missing_index = which(is.na(data.tmp3),arr.ind = T)
  if(length(missing_index) !=0 ){
    for(j in 1:nrow(missing_index)){
      r = missing_index[j,1]
      c = missing_index[j,2]
      data.tmp3[r,c] = mean(c(data.tmp3[r+1,c], data.tmp3[r-1,c], data.tmp3[r,c-1], data.tmp3[r,c+1]), na.rm = T) 
    }
  }
  data.tmp2$dose.response.mats[[1]] = data.tmp3
  
  hsa.tmp = CalculateSynergy2(data.tmp2, method = "HSA", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
  bliss.tmp = CalculateSynergy2(data.tmp2, method = "Bliss", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
  zip.tmp = CalculateSynergy2(data.tmp2, method = "ZIP", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
  loewe.tmp = CalculateSynergy2(data.tmp2, method = "Loewe", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)

  data.tmp$synergy_zip = apply(data.tmp[,c("row","col")],1,function(x) zip.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_hsa = apply(data.tmp[,c("row","col")],1,function(x) hsa.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_bliss = apply(data.tmp[,c("row","col")],1,function(x) bliss.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_loewe = apply(data.tmp[,c("row","col")],1,function(x) loewe.tmp$scores[[1]][x[1],x[2]])
  
  scores[[i]] = data.tmp
  # if (unique(data.tmp$block_id)!=i) 
  #   print(i)
  # flush.console()
}
options(show.error.messages = TRUE)

response_with_scores = do.call(rbind, scores)
write.csv(response_with_scores, "response_oneil_final.csv", row.names = F)
# response_v1 = read.xlsx("response_oneil.xlsx")
# response_v2 = read.xlsx("response_oneil_v2.xlsx")
# cor.test(response_with_scores$synergy_zip, response_v2$synergy_zip)
# cor(response_with_scores[,c("synergy_zip","synergy_bliss","synergy_loewe","synergy_hsa")])
save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Oneil/oneil_part1.RData")

# ---------------------------------------------------
# ask Joseph to provide the summary data
# ---------------------------------------------------


# Now produce the surface.xlsx and curve.xlsx
library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
library(nleqslv)
library(gplots)
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Oneil")
rm(list=ls(all=TRUE))

source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R')

data = read.csv("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Oneil/response_oneil_final.csv", stringsAsFactors = F)
# add row and column numbers
# data2 = ddply(data, c("cell_line_name","drug_row","drug_col","block_id"), transform, row=own_rank(conc_r), col=own_rank(conc_c))

# response as inhibition
data_response <- ReshapeData2(data, data.type = "inhibition")

combo_response = list()
single_response = list()
options(show.error.messages = F)

for(i in unique(data$block_id)){
  cat('\r',i)
  flush.console()
  index = which(data_response$drug.pairs$blockIDs == i)
  
  # surface krigging
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_response$dose.response.mats[[index]])
      # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
      # heatmap.2(mat1,Rowv = F, Colv=F, dendrogram = 'none',trace='none',col = my_palette, density.info="none")
    }, 
    error = function(cond){
      print(i);
      mat1 = data_response$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1)
      return(mat1)
    }
  )
  
  # list format
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
    { 
      tmp = drm(single.row$response ~ single.row$conc_r, fct = LL.4())
      # y = (tmp$coefficients[3] + tmp$coefficients[2] * (conc/tmp$coefficients[4])^tmp$coefficients[1])/(1 + (x/dtmp$coefficients[4])^tmp$coefficients[1])
    }, 
    error = function(cond){
      tmp = drm(single.row$response ~ single.row$conc_r, fct = L.4())
      # y = tmp$coefficients[2] + (tmp$coefficients[3]-tmp$coefficients[2])/(1+exp(tmp$coefficients[1]*(log(x)-log(tmp$coefficients[4]))))
    }
  )
  
  coef.col = tryCatch (
    { 
      tmp = drm(single.col$response ~ single.col$conc_c, fct = LL.4())
    }, 
    error = function(cond){
      tmp = drm(single.col$response ~ single.col$conc_c, fct = L.4())
      # y = tmp$coefficients[2] + (tmp$coefficients[3]-tmp$coefficients[2])/(1+exp(tmp$coefficients[1]*(log(x)-log(tmp$coefficients[4]))))
    }
  )
  
  mat6 = data.frame(matrix(NA, nrow = 2, ncol = 10))
  mat6 = setNames(mat6, c("block_id","drug_row","drug_col","conc_r_unit","conc_c_unit","b","c","d","e","model"))
  
  mat6[1,] =  c(data_response$drug.pairs$blockIDs[index], data_response$drug.pairs$drug.row[index], NA,
                data_response$drug.pairs$concRUnit[index], NA, coef.row$coefficients, as.character(coef.row$call$fct))
  mat6[2,] =  c(data_response$drug.pairs$blockIDs[index], NA, data_response$drug.pairs$drug.col[index], NA,
                data_response$drug.pairs$concCUnit[index], coef.col$coefficients, as.character(coef.col$call$fct))  
  single_response[[i]] = mat6
  
  # fitted.row = coef.row[2] + (coef.row[3]-coef.row[2])/(1+exp(coef.row[1]*(log(x.row)-log(coef.row[4]))))
  # plot(curve.row)
  # points(x.row, fitted.row, col = "red")
}
options(show.error.messages = T)

res_table_response = do.call(rbind, combo_response)
curve = do.call(rbind, single_response)
write.csv(curve, "curve_oneil.csv", row.names = F)

# zip score
colnames(data)[5] = "NA"
colnames(data)[10] = "response"
data_zip<- ReshapeData2(data, data.type = "inhibition")

res_zip = list()
for(i in unique(data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_zip$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_zip$dose.response.mats[[index]])
      # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
      # heatmap.2(mat1,Rowv = F, Colv=F, dendrogram = 'none',trace='none',col = my_palette, density.info="none")
    }, 
    error = function(cond){
      print(i);
      mat1 = data_zip$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_zip'))
  mat2$block_id = data_zip$drug.pairs$blockIDs[i]
  res_zip[[i]] = mat2
}
res_table_zip = do.call(rbind, res_zip)

# bliss score
colnames(data)[10] = "NA"
colnames(data)[11] = "response"
data_bliss<- ReshapeData2(data, data.type = "inhibition")

res_bliss = list()
for(i in unique(data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_bliss$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_bliss$dose.response.mats[[index]])
    }, 
    error = function(cond){
      print(i);
      mat1 = data_bliss$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_bliss'))
  mat2$block_id = data_bliss$drug.pairs$blockIDs[i]
  res_bliss[[i]] = mat2
}
res_table_bliss = do.call(rbind, res_bliss)


# hsa score
colnames(data)[11] = "NA"
colnames(data)[13] = "response"
data_hsa<- ReshapeData2(data, data.type = "inhibition")

res_hsa = list()
for(i in unique(data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_hsa$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_hsa$dose.response.mats[[index]])
    }, 
    error = function(cond){
      print(i);
      mat1 = data_hsa$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_hsa'))
  mat2$block_id = data_hsa$drug.pairs$blockIDs[i]
  res_hsa[[i]] = mat2
}
res_table_hsa = do.call(rbind, res_hsa)


# loewe score
colnames(data)[13] = "NA"
colnames(data)[12] = "response"
data_loewe<- ReshapeData2(data, data.type = "inhibition")

res_loewe = list()
for(i in unique(data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_loewe$drug.pairs$blockIDs == i)
  mat1 = tryCatch (
    { # in case error happens, return the score of NA
      mat1 = smoothing(data_loewe$dose.response.mats[[index]])
    }, 
    error = function(cond){
      print(i);
      mat1 = data_loewe$dose.response.mats[[index]]
      mat1[] = 0
      mat1 = smoothing(mat1)
      return(mat1)
    }
  )
  
  mat2 = setNames(melt(mat1), c('conc_r', 'conc_c', 'synergy_loewe'))
  res_loewe[[i]] = mat2
}
options(show.error.messages = T)
res_table_loewe = do.call(rbind, res_loewe)

surface = cbind(res_table_response, res_table_zip, res_table_bliss, res_table_hsa, res_table_loewe)
write.csv(surface[, c(1:8,11,15,19,23)], "surface_oneil.csv", row.names = F)

save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Oneil/oneil_part2.RData")


