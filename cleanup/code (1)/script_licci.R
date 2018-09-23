library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
library(nleqslv)
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Licci")
rm(list=ls(all=TRUE))

data3 = read.xlsx("nchembio.2382-S3_curated.xlsx")
# POC_Norm_noOutliers is the original viability
# POC_Norm_noOutliers_curated is the harmonized viability

data4 = read.xlsx("nchembio.2382-S4.xlsx")
names(data3) # data3 is all the comb response, while single drug response is missing
names(data4) # data4 is selected drug combo, while single drug response is included 

data3$id = paste(data3$CLOUD_ID_A, data3$CLOUD_ID_B, sep="_")
data4$id = paste(data4$CLOUD.ID1, data4$CLOUD.ID2, sep="_")
data4$POC_Norm_noOutliers = data3$POC_Norm_noOutliers[match(data4$id, data3$id)] 

# confirm that in data4 it is inhibition while in data3 is viability
cor(data4$yc, data4$POC_Norm_noOutliers, use = "complete.obs") # -0.89; not identical but close; how to infer the formula?


# average viablity for each drug in data 3
res1 = ddply(data3, .(CLOUD_ID_A), summarise, mean = mean(POC_Norm_noOutliers))
res2 = ddply(data3, .(CLOUD_ID_B), summarise, mean = mean(POC_Norm_noOutliers))
identical(res1$CLOUD_ID_A,res2$CLOUD_ID_B) # TRUE

cor(res1$mean, res2$mean) # 0.96; safe to take the average viability in drug comb as the estimate of single drug

# compare to inhibition in data4
res3 = ddply(data4, .(CLOUD.ID1), summarise, mean = mean(y1)) # drug1 single inhibition
res4 = ddply(data4, .(CLOUD.ID2), summarise, mean = mean(y2)) # drug2 single inhibition
cor(res3$mean,res4$mean) # 0.99, not identical but close
plot(res3$mean,res3$mean)

res3$mean2 = res1$mean[match(res3$CLOUD.ID1, res1$CLOUD_ID_A)] # viability from data3
cor(res3$mean, res3$mean2, use = "complete.obs") # -0.98
mean(abs(res3$mean-(100-res3$mean2)), na.rm = T) # 11.18

res4$mean2 = res1$mean[match(res4$CLOUD.ID2, res1$CLOUD_ID_A)] # viability from data3
cor(res4$mean, res4$mean2, use = "complete.obs") # -0.97
mean(abs(res4$mean-(100-res4$mean2)), na.rm = T) # 11.29


# so it is safe to use average viability to estimate the single drug sensitivity, from data3
single.viability = cbind(res1, res2$mean, avg = (res1$mean+res2$mean)/2)
plot(data4$yc, data4$POC_Norm_noOutliers)

# now fit the linear model to predict inhibition from viability, based on data4

# linear regression
# fit = lm(yc~POC_Norm_noOutliers, data4)

# loess regression
missing = which(is.na(data4$POC_Norm_noOutliers)==T)
data4.fit = data4[-missing, ]

# fit = loess(yc~POC_Norm_noOutliers, data4.fit) # loess is better
# plot(yc~POC_Norm_noOutliers, data4.fit, pch=19, cex=0.1)
# j <- order(data4.fit$POC_Norm_noOutliers)
# lines(data4.fit$POC_Norm_noOutliers[j], fit$fitted[j], col="red", lwd=3)

# logistic curve, this is the selected model
fit <- drm(yc ~ POC_Norm_noOutliers, data = data4.fit, fct = L.4())
plot(fit)

data_all = list()
for(i in 1:nrow(data3)){
  cat('\r',i)
  drug_row = data3$DrugName_A[i]
  drug_col = data3$DrugName_B[i]
  drug_row_id = data3$CLOUD_ID_A[i]
  drug_col_id = data3$CLOUD_ID_B[i]
  
  conc_r = data3$`SCREENING_A.(µM)`[i]
  conc_c = data3$`SCREENING_B.(µM)`[i]
  conc_r_unit = "uM"
  conc_c_unit = "uM"
  viability = data3$POC_Norm_noOutliers[i]
  
  tmp = expand.grid(block_id = i, conc_r=c(0,conc_r), conc_c=c(0,conc_c))

  index1 = which(tmp$conc_r==0 & tmp$conc_c==0)
  index2 = which(tmp$conc_r==conc_r & tmp$conc_c==0)
  index3 = which(tmp$conc_r==0 & tmp$conc_c==conc_c)
  index4 = which(tmp$conc_r==conc_r & tmp$conc_c==conc_c)
  
  # option1: save the predicted inhibition according to regression model
  tmp$response[index1] = 0
  tmp$response[index2] =
       predict(fit, newdata = data.frame(POC_Norm_noOutliers = single.viability$avg[which(single.viability$CLOUD_ID_A==drug_row_id)]))
  tmp$response[index3] =
       predict(fit, newdata = data.frame(POC_Norm_noOutliers = single.viability$avg[which(single.viability$CLOUD_ID_A==drug_col_id)]))
  tmp$response[index4] = predict(fit, newdata = data.frame(POC_Norm_noOutliers = viability))

  
  # # option2: save the predicted inhibition according to 100-viability
  # tmp$response[index1] = 0
  # tmp$response[index2] = 100 - single.viability$avg[which(single.viability$CLOUD_ID_A==drug_row_id)]
  # tmp$response[index3] = 100 - single.viability$avg[which(single.viability$CLOUD_ID_A==drug_col_id)]
  # tmp$response[index4] = 100 - viability
  
  tmp$drug_row = drug_row
  tmp$drug_col = drug_col
  tmp$conc_r_unit = conc_r_unit
  tmp$conc_c_unit = conc_c_unit
  tmp$drug_row_id = drug_row_id
  tmp$drug_col_id = drug_col_id
  
  data_all[[i]] = tmp
  flush.console()
}
data_licci = do.call(rbind, data_all)
names(data_licci)



# create the template format
template = read.xlsx("template_Licci.xlsx","response",rows = 1)
response_licci = data.frame(matrix(NA, nrow = nrow(data_licci), ncol = ncol(template)))
response_licci = setNames(response_licci, colnames(template))

response_licci$row_id = seq(1:nrow(data_licci))
response_licci$block_id = data_licci$block_id
response_licci$conc_r = data_licci$conc_r
response_licci$conc_c = data_licci$conc_c
response_licci$response = data_licci$response
response_licci$drug_row = data_licci$drug_row
response_licci$drug_col = data_licci$drug_col
response_licci$conc_r_unit = data_licci$conc_r_unit
response_licci$conc_c_unit = data_licci$conc_c_unit
response_licci$drug_row_id = data_licci$drug_row_id
response_licci$drug_col_id = data_licci$drug_col_id
response_licci$cell_line_name = "KBM-7"
response_licci$cell_line_id = "c091"

# QUALITY check
y.data3 = list()
y.data4 = list()
for(i in 1:nrow(data4)){
  test = which(response_licci$drug_row_id==data4$CLOUD.ID1[i] & response_licci$drug_col_id==data4$CLOUD.ID2[i])
  if (length(test)==0) y.data3[[i]] = c(NA, NA, NA) 
  else {
    tmp = response_licci[test, ]
    y1.new = mean(tmp$response[which(tmp$conc_c==0 & tmp$conc_r!=0)]) # avoid CLOFARABINE replicates
    y2.new = mean(tmp$response[which(tmp$conc_c!=0 & tmp$conc_r==0)])
    yc.new = mean(tmp$response[which(tmp$conc_c!=0 & tmp$conc_r!=0)])
    y.data3[[i]] = c(y1.new,y2.new,yc.new)
  }
  y.data4[[i]] = c(data4$y1[i],data4$y2[i],data4$yc[i])
}

y.data3 = do.call(rbind, y.data3)
y.data4 = do.call(rbind, y.data4)

cor(y.data3[,1],y.data4[,1],use="complete.obs") # drug_row consistency option1/2 0.96/0.96
cor(y.data3[,2],y.data4[,2],use="complete.obs") # drug_col consistency 0.95/0.95
cor(y.data3[,3],y.data4[,3],use="complete.obs") # drug comb consistency 0.92/0.89
mean(abs(y.data3[,1]-y.data4[,1]),na.rm=T) # 11.8/12.6
mean(abs(y.data3[,2]-y.data4[,2]),na.rm=T) # 12.3/13.1
mean(abs(y.data3[,3]-y.data4[,3]),na.rm=T) # 9.0/13.8


test = which(response_licci$drug_row_id=="CLOUD117" & response_licci$drug_col_id=="CLOUD171") # the 1st combo in data4
response_licci[test,] # seems fine, the difference is expected.


# ---------------------------------------------------------------
# ------------------- curation stage, not run -------------------
# ---------------------------------------------------------------
# Drug Table
drug_licci = unique(c(response_licci$drug_row, response_licci$drug_col)) # 284 drugs
drug_master = read.csv("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/drug/sucess.csv", stringsAsFactors = F)

# case insensitive search
tmp = lapply(tolower(drug_licci), function(x) grep(x, tolower(drug_master$used.name), fixed = TRUE))
drug_licci[which(sapply(tmp, function(x) length(x) == 0)==T)] # 4 drugs cannot be found
# [1] "Pyridostigmine  Bromide"                 "Teriflunomide (A77 1726)"                "Dideoxycytidine 5'-Triphosphate (Ddctp)"
# [4] "Clofarabine active form" 

# Clofarabine active form (CLOUD282) is merged into Clofarabine (CLOULD209), as CLOUD282 IS NOT FOUND IN nchembio.2382-S2.xlsx
# go back to nchembio.2382-S3.xlsx to change the drug names/done
# rerun the whole analysis/done
# 283 unique drugs
# ----------------- end of curation stage, not run -----------------


# ---------------- Add the Licci drugs to the drug table ------------
drug_licci = unique(c(response_licci$drug_row, response_licci$drug_col)) # 283 drugs
drug_master = read.csv("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/drug/sucess.csv", stringsAsFactors = F)
tmp = lapply(tolower(drug_licci), function(x) which(tolower(drug_master$used.name)==x)) # exact match
drug_licci[which(sapply(tmp, function(x) length(x) == 0)==T)] # 0 drugs cannot be found

drug_table_licci = drug_master[sapply(tmp, function(x) x[1]),c("used.name","Synonym","inchikey","smiles","CID","iupac","molecular_formula","clinical.phase")]
drug_table_licci$used.name = as.character(drug_table_licci$used.name)
colnames(drug_table_licci)[1] = "drug_name"
colnames(drug_table_licci)[2] = "synonym"
colnames(drug_table_licci)[5] = "cid"
colnames(drug_table_licci)[8] = "clinical_phase"

# compare to the existing drug table in the database
drug_table_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\drug.xlsx")
common = intersect(drug_table_licci$cid, drug_table_database$cid) # drugs were found in existing database

tmp = cbind(drug_table_licci$drug_name[match(common, drug_table_licci$cid)], drug_table_database$drug_name[match(common, drug_table_database$cid)])
colnames(tmp) = c("licci","database")
tmp = data.frame(tmp)
identical(tmp$licci, tmp$database) # names should match

drug_table_licci$drug_id = drug_table_database$drug_id[match(drug_table_licci$cid, drug_table_database$cid)]
drug_table_licci$chembl_id = NA
drug_table_licci = drug_table_licci[,c("drug_name","drug_id","synonym","chembl_id","inchikey","smiles","cid","iupac","molecular_formula","clinical_phase")]
write.xlsx(drug_table_licci, "drug_table_licci.xlsx")

# -------------------------------------------------------------
# Quality check
# each row is unique
length(unique(drug_table_licci$cid)) == nrow(drug_table_licci)

# update the drug table
# manual remove []
# copy new drugs into database

# now update the drug ids in the response data
drug_table_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\drug.xlsx")
response_licci$drug_row_id = drug_table_database$drug_id[match(response_licci$drug_row, drug_table_database$drug_name)]
response_licci$drug_col_id = drug_table_database$drug_id[match(response_licci$drug_col, drug_table_database$drug_name)]

# quality check
# PHENPROCOUMON + FLUTAMIDE is the positive findings of the paper
response_licci[which(response_licci$drug_row=="PHENPROCOUMON" & response_licci$drug_col=="FLUTAMIDE"),]

# Save the results
write.xlsx(response_licci,"response_licci_raw.xlsx")


# ----------------------------------------------------------------
# Calculate the synergy scores, use L4 instead of LL4
# Make sure that the data is solid
# ----------------------------------------------------------------
library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
library(nleqslv)

setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Licci")
rm(list=ls(all=TRUE))
response_licci = read.xlsx("response_licci_raw.xlsx")
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/FittingSingleDrug2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/CalculateSynergy2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ZIP2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Loewe2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Bliss2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/HSA2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/BaselineCorrectionSD2.R')

# add row and column numbers
response_licci = ddply(response_licci,c("cell_line_name","drug_row","drug_col","block_id"),transform, row=own_rank(conc_r), col=own_rank(conc_c))

scores = list()
m = length(unique(response_licci$block_id))
options(show.error.messages = F)

for(i in 1:m){
  cat('\r',i)
  index = which(response_licci$block_id == i)
  data.tmp = response_licci[index,]
  data.tmp2 = ReshapeData2(data.tmp, data.type = "inhibition")
  
  hsa.tmp = CalculateSynergy2(data.tmp2, method = "HSA", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
  bliss.tmp = CalculateSynergy2(data.tmp2, method = "Bliss", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
  zip.tmp = CalculateSynergy2(data.tmp2, method = "ZIP", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
  loewe.tmp = CalculateSynergy2(data.tmp2, method = "Loewe", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)

  data.tmp$synergy_zip = apply(data.tmp[,c("row","col")],1,function(x) zip.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_hsa = apply(data.tmp[,c("row","col")],1,function(x) hsa.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_bliss = apply(data.tmp[,c("row","col")],1,function(x) bliss.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_loewe = apply(data.tmp[,c("row","col")],1,function(x) loewe.tmp$scores[[1]][x[1],x[2]])
  
  scores[[i]] = data.tmp
  flush.console()
}
options(show.error.messages = TRUE)

response_licci_with_scores = do.call(rbind, scores)
response_licci_with_scores$block_id = response_licci_with_scores$block_id + 22737
write.csv(response_licci_with_scores,"response_licci_final.csv", row.names = F)
save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Licci/licci_part1.RData")

# ---------------------------------------------------
# ask Joseph to provide the summary.csv data
# ---------------------------------------------------

# ------------------------------------------------
# Now produce the surface.xlsx and curve.xlsx
# ------------------------------------------------
library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Licci")
rm(list=ls(all=TRUE))

source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R')

data = read.csv("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Licci//response_licci_final.csv", stringsAsFactors = F)
# add row and column numbers
# data2 = ddply(data,c("cell_line_name","drug_row","drug_col","block_id"),transform, row=own_rank(conc_r), col=own_rank(conc_c))

# names(data2)[2] = "BlockID"
# names(data2)[3] = "ConcRow"
# names(data2)[4] = "ConcCol"
# names(data2)[5] = "Response"
# names(data2)[6] = "DrugRow"
# names(data2)[7] = "DrugCol"
# names(data2)[8] = "ConcRowUnit"
# names(data2)[9] = "ConcColUnit"

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
  
  # there are two data points, no point to estimate the curve
  # but if insist, use L.4 not LL.4 to avoid errors
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
# options(show.error.messages = T)

res_table_response = do.call(rbind, combo_response)
curve = do.call(rbind, single_response)

curve$block_id = as.numeric(curve$block_id)
write.csv(curve, "curve_licci.csv", row.names = F)

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
res_table_loewe = do.call(rbind, res_loewe)

surface = cbind(res_table_response, res_table_zip, res_table_bliss, res_table_hsa, res_table_loewe)

# surface$block_id = surface$block_id
write.csv(surface[,c(1:8,11,15,19,23)], "surface_licci.csv", row.names = F)

save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Licci/licci_part2.RData")


