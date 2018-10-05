# MIT license
# this file uses NCI ALMANAC data as input: 
# https://wiki.nci.nih.gov/download/attachments/338237347/ComboDrugGrowth_Nov2017.zip?version=1&modificationDate=1510057275000&api=v2
# it reshapes and modifies to be used fpr synergyfinder::ReshapeData()
# authors (alphabetical order): Saad, J., Tang, J., Zagidullin, B. 
# version 0.51

#working dir laptop
setwd('/Users/zagidull/Documents/git/sharing_notebooks/cleanup')
#working dir server
#setwd('/Users/zagidull/Desktop/synergy_calc_august/cleanup')


#housekeeping
rm(list=ls(all=TRUE))

#required libraries
library("dplyr")
library('drc')
library('nleqslv')
library('readr')
library('reshape2')
library('tibble')

#required functions
source('correctingConc.R')
source('rowInsert.R')
source('ReshapeData.R')
source('exchange.R')
source('ZIP.R')
source('HSA.R')
source('Bliss.R')
source('Loewe.R')
source('BaselineCorrectionSD2.R')
source('FittingSingleDrug.R')
source('CalculateSynergy.R')
source('ReshapeForDB.R')
source('ggplotRegression.R')


#input file
file <- c('ComboDrugGrowth_Nov2017.csv')

#reading in the data and correcting col classes
input.data <- read_csv(file, progress = T) # it is nice to read as tibble, this helps get rid of formatting errors

input.data$TESTDATE <- as.Date(input.data$TESTDATE,"%m/%d/%Y") #correcting DATA col class
#sapply(temp.for_sorting, class) this checks the class of each of the columns
input.data -> temp.for_sorting

# introduce a new column "NATREATMENT". it is a factor. With levels (noTZ_to_noTZ, mean_to_noTZ, TZ_to_noTZ)
#temp.for_sorting$NATREATMENT <- factor(x = rep(x = NA, times = nrow(temp.for_sorting)),
#                                       levels = c('noTZ_to_noTZ', 'mean_to_noTZ', 'TZ_to_noTZ'), 
#                                       ordered = T)

# create a numeric column, where 1 is noTZ_to_noTZ, 2 is TZ_to_noTZ, 3 is mean_to_noTZ
# temp.for_sorting$NATREATMENT <- as.numeric(NA)

#creating placeholde for results  
reshaped.input.datalist <- list()

#cellnames to be searched for
temp.cellnames_old <- unique(temp.for_sorting$CELLNAME)
temp.cellnames <- temp.cellnames_old[-61] #todo check why the last cell line has to be removed

# random sampling of cell lines
# set.seed(1)
# temp.cellnames <- sample(temp.cellnames, 3, replace = F)

# there are 208605 rows where PERCENTINHIBITIONNOTZ is NA
# Substitute PERCENTGROWTHNOTZ with PERCENTGROWTH, when PERCENTGROWTHNOTZ == NA
# https://stackoverflow.com/questions/4862178/remove-rows-with-all-or-some-nas-missing-values-in-data-frame
#temp.for_sorting[complete.cases(temp.for_sorting[, 21]), ] -> temp.for_sorting

# all the rows where PERCENTGROWTHNOTZ is NOT na. Which means we do nothing and set factor as noTZ_to_noTZ
# temp.for_sorting[complete.cases(temp.for_sorting[, 21]), ]$NATREATMENT <- 1

# all the rows where there is NA in PERCENTGRWOTHNOTZ and PERCENTGROWTH >= 0
# temp.for_sorting[ (is.na(temp.for_sorting$PERCENTGROWTHNOTZ) & temp.for_sorting$PERCENTGROWTH >= 0), ]$NATREATMENT <- 2

# all the rows where there is NA in PERCENTGRWOTHNOTZ and PERCENTGROWTH < 0
# temp.for_sorting[ (is.na(temp.for_sorting$PERCENTGROWTHNOTZ) & temp.for_sorting$PERCENTGROWTH < 0), ]$NATREATMENT <- 3

# fit linear model. PERCENTGROWTHNOTZ - response, PERCENTGROWTH - predictor. Trying fitting logistic regression, so x is from -inf to +inf and y is always positive
temp.noTZ_from_TZ <- lm(formula = PERCENTGROWTHNOTZ ~ PERCENTGROWTH, data = temp.for_sorting, na.action=na.omit)
saveRDS(file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "lm_model"), object = temp.noTZ_from_TZ)

#intercept <- 35.0813
#fit1 <- lm(I(PERCENTGROWTHNOTZ - intercept) ~ 0 + PERCENTGROWTH, data = temp.for_sorting, na.action=na.omit)
#fit2 <- lm(formula = PERCENTGROWTHNOTZ ~ 0 + PERCENTGROWTH, offset = rep(intercept, length(PERCENTGROWTHNOTZ)), data = temp.for_sorting, na.action=na.omit)

# for graphing
# d <- ggplot(fit2$model, aes(x = PERCENTGROWTH, y = PERCENTGROWTHNOTZ))
# d + stat_binhex()
# ggplotRegression(noTZ_from_TZ) + geom_abline(intercept = 35.08, slope = 0.6333)


# predicting for when PERCENTGROWTHNOTZ is NA
temp.to.predict <- temp.for_sorting[ (is.na(temp.for_sorting$PERCENTGROWTHNOTZ)),]
temp.predicted <- stats::predict(temp.noTZ_from_TZ, temp.to.predict)

# let's set PERCENTGROWTHNOTZ < 0, as 0. The significance is that it works as good as the Positive Control. I guess
temp.predicted[temp.predicted < 0] <- 0

#exchanging NA in PERCENTGROWTHNOTZ with predicted ones
temp.for_sorting[ (is.na(temp.for_sorting$PERCENTGROWTHNOTZ)),]$PERCENTGROWTHNOTZ <- temp.predicted

#these are still negative PERECENTGROWTHNOTZ
temp.for_sorting[temp.for_sorting$PERCENTGROWTHNOTZ < 0,] -> temp.negative #9080 rows. Set them as zero
temp.for_sorting[temp.for_sorting$PERCENTGROWTHNOTZ >= 0,] -> temp.positive #3677395 rows

if ((nrow(temp.positive) == nrow(input.data)) && nrow(temp.negative) == 0 ) {
  cat('\r','NB: no more negative values of percent growth are present.Size of processed file is identical to input file')
  flush.console()
}

for (temp.name in temp.cellnames) {
  #this creates a place holder for one cellline 
  temp.for_sorting_many <- list()
  temp.for_sorting_many <- temp.for_sorting %>% # taking one cell line as an example////nb this is the original input file
    filter(CELLNAME==temp.name)
  
  #this makes a list of drugs used in given cellline
  temp.uniqueDrug1Drub2Combos <- list()
  temp.uniqueDrug1Drub2Combos <- temp.for_sorting_many %>%
    mutate(NSC1,NSC2) %>%
    dplyr::select(NSC1, NSC2) %>%
    unique %>%
    na.omit()
  
  # for now checking with a randomly selected set of drug combos
  # set.seed(1)
  # temp.uniqueDrug1Drub2Combos <- temp.uniqueDrug1Drub2Combos[sample(nrow(temp.uniqueDrug1Drub2Combos), 5),]
  
  for (temp.i in 1:nrow(temp.uniqueDrug1Drub2Combos)) {
    temp.perrow <- list()
    temp.input <- list()
    temp.output <- list()
    temp.param1 <- vector(mode="numeric", length=1)
    temp.param2 <- vector(mode="numeric", length=1)
    temp.perrow<- temp.uniqueDrug1Drub2Combos[temp.i,]
    
    temp.param1 <- as.numeric(temp.perrow[1])
    temp.param2 <- as.numeric(temp.perrow[2])
    
    temp.input <- correctingConc(temp.for_sorting_many, temp.param1, temp.param2)
    temp.iden <- paste(temp.i, temp.name, sep = ":")
    temp.output <-  tryCatch(rowInsert(temp.input,temp.iden), error = function(e) NA)
    #    output1 <- rowInsert(input,i)
    reshaped.input.datalist[[temp.name]][[temp.i]] <- temp.output
    if (temp.i%%1000==0) cat('left', nrow(temp.uniqueDrug1Drub2Combos) - temp.i,'\n',sep='' )
    if (temp.i == nrow(temp.uniqueDrug1Drub2Combos)) cat('cell group processed = ', temp.name, ' \n', sep='')
    
  }
  
}

# bind all in one and leave only rows with the correct number of columns. So 3x3 + 3 + 3 + 1 = 16 or 3x5 + 3 + 5 + 1 = 24
do.call(rbind, do.call(rbind, reshaped.input.datalist)) -> temp.unbound
temp.unbound.filtered <- temp.unbound %>%
  dplyr::group_by(BlockID) %>%
    dplyr::mutate(nrows = n()) %>%
      dplyr::filter(nrows == 16 | nrows == 24)

# these were filtered out? 
# dplyr::setdiff(temp.unbound,temp.unbound.filtered[,-length(temp.unbound.filtered)]) -> testing
# unique(lapply(testing$BlockID, function (x) unlist(strsplit(x, ":"))[2])) -> to check which CELLLINES are used




#renames colnames
##tracemem says that by renaming columns we make it a completely new object. Which means copying over
##tracemem(unbound_filtered) -> before
names(temp.unbound.filtered) -> temp.saving.names
names(temp.unbound.filtered) <- c('row', 'col', 'block_id', 'response', 'conc_c_unit', 'conc_r_unit', 'conc_r', 'conc_c', 'drug_row', 'drug_col', 'nrows')
##tracemem(unbound_filtered) -> after

#reshaping
temp.reshaped <- ReshapeData(temp.unbound.filtered, data.type = 'viability')

#change concetrations from M to uM
temp.reshaped$dose.response.mats <- lapply(temp.reshaped$dose.response.mats, function(x) {
  colnames(x) <- as.numeric(colnames(x))*(10^6)
  rownames(x) <- as.numeric(rownames(x))*(10^6)
  return(x)
})
temp.reshaped$drug.pairs$concRUnit <- temp.reshaped$drug.pairs$concCUnit <- 'uM'

# the idea is to wrap CalculateSynergy and do the reshape for DB as a function. Ideally, as a function that could be parallelized.
# could be either screen or via parallel R package
#temp.methods <- c('ZIP','Loewe','Bliss','HSA')
#for (temp.method in temp.methods){
#  CalculateSynergy(reshaped, method = temp.method, correction = T, Emin = 0, Emax = NA)
#}

#saving to run on server
saveRDS(object = temp.reshaped, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped"))

#used for testing. getting first 10 drug pairs per each cell line
#temp.reshaped -> temp.old.reshaped
#temp.reshaped$dose.response.mats[c(1:10)] -> temp.reshaped$dose.response.mats
#temp.reshaped$drug.pairs[c(1:10),] -> temp.reshaped$drug.pairs

# calculating synergy. reshaped.* files are saved on disk
CalculateSynergy(temp.reshaped, method = 'ZIP', correction = T, Emin = 0, Emax = NA) -> temp.reshaped.ZIP
CalculateSynergy(temp.reshaped, method = 'Loewe', correction = T, Emin = 0, Emax = NA) -> temp.reshaped.Loewe
CalculateSynergy(temp.reshaped, method = 'Bliss', correction = T, Emin = 0) -> temp.reshaped.Bliss
CalculateSynergy(temp.reshaped, method = 'HSA', correction = T, Emin = 0, Emax = NA) -> temp.reshaped.HSA

# saving after reshape
saveRDS(object = temp.reshaped.ZIP, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_ZIP"))
saveRDS(object = temp.reshaped.Loewe, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_Loewe"))
saveRDS(object = temp.reshaped.Bliss, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_Bliss"))
saveRDS(object = temp.reshaped.HSA, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_HSA"))

# reshaping for kriegging and CSS calc
do.call(rbind, ReshapeForDB(temp.reshaped.ZIP, 'ZIP')) -> temp.rbound.ZIP
do.call(rbind, ReshapeForDB(temp.reshaped.Loewe, 'Loewe')) -> temp.rbound.Loewe
do.call(rbind, ReshapeForDB(temp.reshaped.Bliss, 'Bliss')) -> temp.rbound.Bliss
do.call(rbind, ReshapeForDB(temp.reshaped.HSA, 'HSA')) -> temp.rbound.HSA

temp.out <-  merge(merge(temp.rbound.ZIP, temp.rbound.Loewe, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID', 'DrugRow', 'DrugCol','DrugColUnit', 'DrugRowUnit')), 
              merge(temp.rbound.HSA, temp.rbound.Bliss, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID', 'DrugRow', 'DrugCol','DrugColUnit', 'DrugRowUnit')), 
              by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID', 'DrugRow', 'DrugCol','DrugColUnit', 'DrugRowUnit'))

saveRDS(object = temp.out, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_All"))

#housekeeping
# rm(list=ls(pattern='^temp'))

#data <- temp.rbound.HSA
# data <- temp.rbound.HSA[temp.rbound.HSA$BlockID == '1:ACHN',]
names(temp.out) -> temp.names.old
temp.data <- temp.out
names(temp.data) <- c('conc_r','conc_c','response','block_id','drug_row','drug_col','conc_c_unit',"conc_r_unit",'ZIP','Loewe', 'HSA', 'Bliss')

# laptop
setwd('/Users/zagidull/Documents/git/sharing_notebooks/cleanup/code (1)')
# server


library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
library(nleqslv)
library(gplots)

source('smoothing.R')
source('own_rank.R')
source('ReshapeData2.R')

source('FittingSingleDrug2.R')
source('FittingSingleDrug2.R')
source('CalculateSynergy2.R')
source('ZIP2.R')
source('Loewe2.R')
source('Bliss2.R')
source('HSA2.R')
source('BaselineCorrectionSD2.R')

temp.response = ddply(temp.data, c("block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))

temp.data_response <- ReshapeData2(temp.response, data.type = "inhibition")

temp.scores = list()
m = length(unique(temp.response$block_id))
options(show.error.messages = F)

temp.combo_response = list()
temp.single_response = list()

for(i in unique(temp.data$block_id)){
  # 10:ACHN comes before 1:ACHN
  cat('\r',i)
  flush.console()
  index = which(temp.data_response$drug.pairs$blockIDs == i)
  
  # surface krigging
  temp.mat1 = tryCatch (
    { # in case error happens, return the score of NA
      temp.mat1 = smoothing(temp.data_response$dose.response.mats[[index]])
      # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
      # heatmap.2(mat1,Rowv = F, Colv=F, dendrogram = 'none',trace='none',col = my_palette, density.info="none")
    }, 
    error = function(cond){
      print(i);
      temp.mat1 = temp.data_response$dose.response.mats[[index]]
      temp.mat1[] = 0
      temp.mat1 = smoothing(temp.mat1)
      return(temp.mat1)
    }
  )
  
  # list format
  temp.mat2 = setNames(melt(temp.mat1), c('conc_r', 'conc_c', 'response'))
  temp.mat2$block_id = temp.data_response$drug.pairs$blockIDs[index]
  temp.mat2$drug_row = temp.data_response$drug.pairs$drug.row[index]
  temp.mat2$drug_col = temp.data_response$drug.pairs$drug.col[index]
  temp.mat2$conc_r_unit = temp.data_response$drug.pairs$concRUnit[index]
  temp.mat2$conc_c_unit = temp.data_response$drug.pairs$concCUnit[index]
  temp.combo_response[[index]] = temp.mat2
  
  # single drug curve fitting
  temp.mat3 = temp.data_response$dose.response.mats[[index]]
  temp.mat4 = setNames(melt(temp.mat3), c('conc_r', 'conc_c', 'response'))
  
  temp.mat5 = temp.mat4[which(temp.mat4$conc_r==0 | temp.mat4$conc_c==0),]
  single.row = temp.mat5[which(temp.mat5$conc_c==0),]
  single.col = temp.mat5[which(temp.mat5$conc_r==0),]
  
  coef.row = tryCatch (
    { 
      tmp = drm(single.row$response ~ single.row$conc_r, fct = LL.4())
      # y1 = (tmp$coefficients[3] + tmp$coefficients[2] * (x/tmp$coefficients[4])^tmp$coefficients[1])/(1 + (x/tmp$coefficients[4])^tmp$coefficients[1])
      # y2 = tmp$coefficients[2] + (tmp$coefficients[3] - tmp$coefficients[2])/(1+exp(tmp$coefficients[1]*(log(x)-log(tmp$coefficients[4]))))
    }, 
    error = function(cond){
      tmp = drm(single.row$response ~ single.row$conc_r, fct = L.4())
      # y1 = tmp$coefficients[2] + (tmp$coefficients[3]-tmp$coefficients[2])/(1+exp(tmp$coefficients[1]*(x-tmp$coefficients[4])))
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
  
  temp.mat6 = data.frame(matrix(NA, nrow = 2, ncol = 10))
  temp.mat6 = setNames(temp.mat6, c("block_id","drug_row","drug_col","conc_r_unit","conc_c_unit","b","c","d","e","model"))
  
  temp.mat6[1,] =  c(temp.data_response$drug.pairs$blockIDs[index], temp.data_response$drug.pairs$drug.row[index], NA,
                temp.data_response$drug.pairs$concRUnit[index], NA, coef.row$coefficients, as.character(coef.row$call$fct))
  temp.mat6[2,] =  c(temp.data_response$drug.pairs$blockIDs[index], NA, temp.data_response$drug.pairs$drug.col[index], NA,
                temp.data_response$drug.pairs$concCUnit[index], coef.col$coefficients, as.character(coef.col$call$fct))  
  temp.single_response[[index]] = temp.mat6
  
  # fitted.row = coef.row[2] + (coef.row[3]-coef.row[2])/(1+exp(coef.row[1]*(log(x.row)-log(coef.row[4]))))
  # plot(curve.row)
  # points(x.row, fitted.row, col = "red")
}
options(show.error.messages = T)

res_table_response = do.call(rbind, temp.combo_response)
curve = do.call(rbind, temp.single_response)

# hsa score
colnames(temp.data)[3] = "NA"
colnames(temp.data)[11] = "response"
temp.response = ddply(temp.data, c("block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))

data_hsa<- ReshapeData2(temp.response, data.type = "inhibition")

res_hsa = list()
for(i in unique(temp.data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_hsa$drug.pairs$blockIDs == i)
  temp.mat1 = tryCatch (
    { # in case error happens, return the score of NA
      temp.mat1 = smoothing(data_hsa$dose.response.mats[[index]])
    }, 
    error = function(cond){
      print(i);
      temp.mat1 = data_hsa$dose.response.mats[[index]]
      temp.mat1[] = 0
      temp.mat1 = smoothing(temp.mat1)
      return(temp.mat1)
    }
  )
  
  temp.mat2 = setNames(melt(temp.mat1), c('conc_r', 'conc_c', 'synergy_hsa'))
  temp.mat2$block_id = data_hsa$drug.pairs$blockIDs[index]
  res_hsa[[index]] = temp.mat2
}
res_table_hsa = do.call(rbind, res_hsa)

# zip score
colnames(temp.data)[11] = "NA"
colnames(temp.data)[9] = "response"
temp.response = ddply(temp.data, c("block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))

data_zip<- ReshapeData2(temp.response, data.type = "inhibition")

res_zip = list()
for(i in unique(temp.data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_zip$drug.pairs$blockIDs == i)
  temp.mat1 = tryCatch (
    { # in case error happens, return the score of NA
      temp.mat1 = smoothing(data_zip$dose.response.mats[[index]])
      # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
      # heatmap.2(temp.mat1,Rowv = F, Colv=F, dendrogram = 'none',trace='none',col = my_palette, density.info="none")
    }, 
    error = function(cond){
      print(i);
      temp.mat1 = data_zip$dose.response.mats[[index]]
      temp.mat1[] = 0
      temp.mat1 = smoothing(temp.mat1)
      return(temp.mat1)
    }
  )
  
  temp.mat2 = setNames(melt(temp.mat1), c('conc_r', 'conc_c', 'synergy_zip'))
  temp.mat2$block_id = data_zip$drug.pairs$blockIDs[index]
  res_zip[[index]] = temp.mat2
}
res_table_zip = do.call(rbind, res_zip)

# bliss score
colnames(temp.data)[9] = "NA"
colnames(temp.data)[12] = "response"
temp.response = ddply(temp.data, c("block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))

data_bliss<- ReshapeData2(temp.response, data.type = "inhibition")

res_bliss = list()
for(i in unique(temp.data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_bliss$drug.pairs$blockIDs == i)
  temp.mat1 = tryCatch (
    { # in case error happens, return the score of NA
      temp.mat1 = smoothing(data_bliss$dose.response.mats[[index]])
    }, 
    error = function(cond){
      print(i);
      temp.mat1 = data_bliss$dose.response.mats[[index]]
      temp.mat1[] = 0
      temp.mat1 = smoothing(temp.mat1)
      return(temp.mat1)
    }
  )
  
  temp.mat2 = setNames(melt(temp.mat1), c('conc_r', 'conc_c', 'synergy_bliss'))
  temp.mat2$block_id = data_bliss$drug.pairs$blockIDs[index]
  res_bliss[[index]] = temp.mat2
}
res_table_bliss = do.call(rbind, res_bliss)

# loewe score
colnames(temp.data)[12] = "NA"
colnames(temp.data)[10] = "response"
temp.response = ddply(temp.data, c("block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))

data_loewe<- ReshapeData2(temp.response, data.type = "inhibition")

res_loewe = list()
for(i in unique(temp.data$block_id)){
  cat('\r',i)
  flush.console()
  
  index = which(data_loewe$drug.pairs$blockIDs == i)
  temp.mat1 = tryCatch (
    { # in case error happens, return the score of NA
      temp.mat1 = smoothing(data_loewe$dose.response.mats[[index]])
    }, 
    error = function(cond){
      print(i);
      temp.mat1 = data_loewe$dose.response.mats[[index]]
      temp.mat1[] = 0
      temp.mat1 = smoothing(temp.mat1)
      return(temp.mat1)
    }
  )
  
  temp.mat2 = setNames(melt(temp.mat1), c('conc_r', 'conc_c', 'synergy_loewe'))
  res_loewe[[index]] = temp.mat2
}
options(show.error.messages = T)
res_table_loewe = do.call(rbind, res_loewe)

surface = cbind(res_table_response, res_table_zip, res_table_bliss, res_table_hsa, res_table_loewe)

rm(list=ls(pattern='^temp'))

