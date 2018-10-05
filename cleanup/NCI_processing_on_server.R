# MIT license
# this file uses NCI ALMANAC data as input: 
# https://wiki.nci.nih.gov/download/attachments/338237347/ComboDrugGrowth_Nov2017.zip?version=1&modificationDate=1510057275000&api=v2
# it reshapes and modifies to be used fpr synergyfinder::ReshapeData()
# authors (alphabetical order): Saad, J., Tang, J., Zagidullin, B. 
# version 0.51

#working dir server
setwd('/home/bulat/NCI/cleanup')


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

CalculateSynergy(temp.reshaped, method = 'HSA', correction = T, Emin = 0, Emax = NA) -> temp.reshaped.HSA
saveRDS(object = temp.reshaped.HSA, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_HSA"))

