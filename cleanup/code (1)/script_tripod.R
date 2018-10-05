library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
library(nleqslv)
# # Preprocessing - done ------------------------------------------
# response_forcina = read.csv("response.csv", stringsAsFactors = F)
# drug_forcina = read.xlsx("drug.xlsx")
# cell_line_forcina = read.xlsx("cell_line.xlsx")
# study_forcina = read.xlsx("study.xlsx")
# assay_forcina = read.xlsx("assay.xlsx")
# 
# drug_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\drug.xlsx")
# cell_line_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\cell_line.xlsx")
# study_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\study.xlsx")
# assay_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\assay.xlsx")
# 
# # QC
# length(intersect(drug_database$drug_id, drug_forcina$drug_id)) == 0 # no overlapped drugs
# length(intersect(drug_database$cid, drug_forcina$cid)) == 0
# 
# # each row is unique
# length(unique(drug_forcina$cid)) == nrow(drug_forcina)
# length(intersect(cell_line_database$cellosaurus_assession, cell_line_forcina$cellosaurus_assession)) == 0 # FALSE The cell line (T98G) is already there
# cell_line_database[match(cell_line_forcina$cellosaurus_assession, cell_line_database$cellosaurus_assession),]
# 
# length(intersect(drug_database$drug_id, drug_forcina$drug_id)) == 0
# length(setdiff(response_forcina$drug_row_id, c(drug_database$drug_id, drug_forcina$drug_id))) == 0
# length(setdiff(response_forcina$drug_col_id, c(drug_database$drug_id, drug_forcina$drug_id))) == 0
# # END QC

# ----------------------------------------------------------------
# Calculate the synergy scores
# Make sure that the data is solid
# ----------------------------------------------------------------
library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
library(nleqslv)
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Griner\\upload_tripod")
rm(list=ls(all=TRUE))

# drug_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\drug.xlsx")
response = read.csv("response_griner.csv", stringsAsFactors = F)
response_shuyu =  read.csv("response_upload.csv", stringsAsFactors = F)

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
response = ddply(response,c("cell_line_name","drug_row","drug_col","block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))

scores = list()
m = unique(response$block_id)
options(show.error.messages = F)

for(i in m){
  cat('\r',i)
  index = which(response$block_id == i)
  data.tmp = response[index,]
  data.tmp2 = ReshapeData2(data.tmp, data.type = "inhibition")
  
  hsa.tmp = CalculateSynergy2(data.tmp2, method = "HSA", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
  bliss.tmp = CalculateSynergy2(data.tmp2, method = "Bliss", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
  zip.tmp = CalculateSynergy2(data.tmp2, method = "ZIP", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
  loewe.tmp = CalculateSynergy2(data.tmp2, method = "LOEWE", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
  
  data.tmp$synergy_zip = apply(data.tmp[,c("row","col")],1,function(x) zip.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_hsa = apply(data.tmp[,c("row","col")],1,function(x) hsa.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_bliss = apply(data.tmp[,c("row","col")],1,function(x) bliss.tmp$scores[[1]][x[1],x[2]])
  data.tmp$synergy_loewe = apply(data.tmp[,c("row","col")],1,function(x) loewe.tmp$scores[[1]][x[1],x[2]])
  
  scores[[i]] = data.tmp
  flush.console()
}
options(show.error.messages = TRUE)

response_with_scores = do.call(rbind, scores)
response_with_scores$block_id = response_with_scores$block_id 
write.csv(response_with_scores,"C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Griner/tang/response_with_scores.csv", row.names = F)
save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/Griner/tang/response_with_scores.RData")

# ---------------------------------------------------
# ask Joseph to provide the summary.csv data
# ---------------------------------------------------
# "summary_forcina.csv"

# ------------------------------------------------
# Now produce the surface.xlsx and curve.xlsx
# ------------------------------------------------
library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\forcina")
rm(list=ls(all=TRUE))

source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R')

data = read.csv("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/forcina//response_forcina.csv", stringsAsFactors = F)
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
options(show.error.messages = T)

res_table_response = do.call(rbind, combo_response)
curve = do.call(rbind, single_response)

curve$block_id = as.numeric(curve$block_id)
write.csv(curve, "curve_forcina.csv", row.names = F)

# zip score
colnames(data)[which(colnames(data) == "response")] = "NA"
colnames(data)[which(colnames(data) == "synergy_zip")] = "response"
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
colnames(data)[which(colnames(data) == "response")] = "NA"
colnames(data)[which(colnames(data) == "synergy_bliss")] = "response"
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
colnames(data)[which(colnames(data) == "response")] = "NA"
colnames(data)[which(colnames(data) == "synergy_hsa")] = "response"
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
colnames(data)[which(colnames(data) == "response")] = "NA"
colnames(data)[which(colnames(data) == "synergy_loewe")] = "response"
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
write.csv(surface[,c(1:8,11,15,19,23)], "surface_forcina.csv", row.names = F)

save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/forcina/forcina_part2.RData")


# # QC before upload to database
# setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\")
# rm(list=ls(all=TRUE))
# response_forcina = read.csv("response_forcina.csv", stringsAsFactors = F)
# drug_forcina = read.xlsx("drug_forcina.xlsx")
# cell_line_forcina = read.xlsx("cell_line_forcina.xlsx")
# study_forcina = read.xlsx("study_forcina.xlsx")
# assay_forcina = read.xlsx("assay_forcina.xlsx")
# curve_forcina = read.csv("curve_forcina.csv", stringsAsFactors = F)
# surface_forcina = read.csv("surface_forcina.csv", stringsAsFactors = F)
# 
# drug_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\drug.xlsx")
# cell_line_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\cell_line.xlsx")
# study_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\study.xlsx")
# assay_database = read.xlsx("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\assay.xlsx")
# 
# # check that the drugs in reponse are either in drug_forcina or drug_database
# length(setdiff(c(unique(response_forcina$drug_row_id), unique(response_forcina$drug_col_id)), c(drug_forcina$id, drug_database$drug_id))) == 0
# # setdiff(c(unique(surface_forcina$drug_row_id), unique(surface_forcina$drug_col_id)), c(drug_forcina$id, drug_database$drug_id))
# length(setdiff(unique(response_forcina$block_id), unique(curve_forcina$block_id))) == 0
# 
# # drug names should match between resonse and drug
# drug_name_response = data.frame(cbind(drug_forcina$name[match(response_forcina$drug_row_id, drug_forcina$id)], drug_database$drug_name[match(response_forcina$drug_row_id, drug_database$drug_id)], 
#                                  toupper(response_forcina$drug_row)), stringsAsFactors = F)
# drug_name_response[is.na(drug_name_response[,1]),1] = drug_name_response[is.na(drug_name_response[,1]),2]
# colnames(drug_name_response) = c("forcina_drug","database","forcina_response")
# identical(drug_name_response[,1], drug_name_response[,3])
# 
# # drug names should match between surface and drug
# drug_name_surface = data.frame(cbind(drug_forcina$name[match(surface_forcina$drug_row, drug_forcina$name)], drug_database$drug_name[match(surface_forcina$drug_row, drug_database$drug_name)], 
#                                  toupper(surface_forcina$drug_row)), stringsAsFactors = F)
# drug_name_surface[is.na(drug_name_surface[,1]),1] = drug_name_surface[is.na(drug_name_surface[,1]),2]
# colnames(drug_name_surface) = c("forcina_drug","database","forcina_surface")
# identical(drug_name_surface[,1], drug_name_surface[,3])


# Now start formating the data tables according to the schema
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\")
rm(list=ls(all=TRUE))
response_forcina = read.csv("response_forcina.csv", stringsAsFactors = F)
drug_forcina = read.xlsx("drug_forcina.xlsx")
cell_line_forcina = read.xlsx("cell_line_forcina.xlsx")
study_forcina = read.xlsx("study_forcina.xlsx")
assay_forcina = read.xlsx("assay_forcina.xlsx")
curve_forcina = read.csv("curve_forcina.csv", stringsAsFactors = F)
surface_forcina = read.csv("surface_forcina.csv", stringsAsFactors = F)
summary_forcina = read.csv("summary_forcina.csv", stringsAsFactors = F)
tissue_forcina = read.xlsx("tissue_forcina.xlsx")
disease_forcina = read.xlsx("disease_forcina.xlsx")

# tissue
write.csv(tissue_forcina, "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\tissue_forcina_upload.csv", row.names = F)
write.csv(disease_forcina, "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\disease_forcina_upload.csv", row.names = F)

# response
head(response_forcina)
write.csv(response_forcina[,c("block_id","conc_r","conc_c","response","synergy_zip","synergy_bliss","synergy_loewe","synergy_hsa")], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\response_forcina_upload.csv", row.names = F)

# cell_line
head(cell_line_forcina)
write.csv(cell_line_forcina, "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\cell_line_forcina_upload.csv", row.names = F)

# drug
names(drug_forcina)
write.csv(drug_forcina, "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\drug_forcina_upload.csv", row.names = F)

# surface
names(surface_forcina)
write.csv(surface_forcina[, -which(names(surface_forcina) %in% c("drug_row","drug_col","conc_r_unit","conc_c_unit"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\surface_forcina_upload.csv", row.names = F)

# summary
names(summary_forcina)
# add drug_row_id and drug_col_id
summary_forcina$drug_row_id = response_forcina$drug_row_id[match(summary_forcina$drug_row, response_forcina$drug_row)]
summary_forcina$drug_col_id = response_forcina$drug_col_id[match(summary_forcina$drug_col, response_forcina$drug_col)]
summary_forcina$cell_line_id = response_forcina$cell_line_id[match(summary_forcina$cell_line_name, response_forcina$cell_line_name)]
write.csv(summary_forcina[, -which(names(summary_forcina) %in% c("drug_row","drug_col","conc_r_unit","conc_c_unit","cell_line_name"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\summary_forcina_upload.csv", row.names = F)

# curve
curve_forcina$drug_row_id = response_forcina$drug_row_id[match(curve_forcina$drug_row, response_forcina$drug_row)]
curve_forcina$drug_col_id = response_forcina$drug_col_id[match(curve_forcina$drug_col, response_forcina$drug_col)]
curve_forcina$drug_row_id[is.na(curve_forcina$drug_row_id)] = 0
curve_forcina$drug_col_id[is.na(curve_forcina$drug_col_id)] = 0
write.csv(curve_forcina[, -which(names(curve_forcina) %in% c("drug_row","drug_col","conc_r_unit","conc_c_unit"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\curve_forcina_upload.csv", row.names = F)

# study
head(study_forcina)
write.csv(study_forcina, "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\study_forcina_upload.csv", row.names = F)

# assay
head(assay_forcina)
colnames(assay_forcina)[which(colnames(assay_forcina) == "assay_type")] = "type"
write.csv(assay_forcina, "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Forcina\\Upload\\assay_forcina_upload.csv", row.names = F)

