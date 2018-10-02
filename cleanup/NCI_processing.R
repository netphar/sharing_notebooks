# MIT license
# this file uses NCI ALMANAC data as input: 
# https://wiki.nci.nih.gov/download/attachments/338237347/ComboDrugGrowth_Nov2017.zip?version=1&modificationDate=1510057275000&api=v2
# it reshapes and modifies to be used fpr synergyfinder::ReshapeData()
# authors: Tang, J., Zagidullin, B.
# version 0.4

# NB: do not forget to modify file

#working dir
setwd('/Users/zagidull/Desktop/synergy_calc_august/cleanup')

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


#input file
file <- c('ComboDrugGrowth_Nov2017.csv')
# new.csv has 9999 lines? so, for like initial testing
#file <- c('new.csv') 

#reading in the data and correcting col classes
input.data <- read_csv(file, progress = F) # it is nice to read as tibble, this helps get rid of formatting errors

input.data$TESTDATE <- as.Date(input.data$TESTDATE,"%m/%d/%Y") #correcting DATA col class
#sapply(temp.for_sorting, class) this checks the class of each of the columns
input.data -> temp.for_sorting

#creating placeholde for results  
reshaped.input.datalist <- list()

#cellnames to be searched for
temp.cellnames_old <- unique(temp.for_sorting$CELLNAME)
temp.cellnames <- temp.cellnames_old[-61] #todo check why the last cell line has to be removed

#remove rows where PERCENTGROWTHNOTZ is na. 
# https://stackoverflow.com/questions/4862178/remove-rows-with-all-or-some-nas-missing-values-in-data-frame
temp.for_sorting[complete.cases(temp.for_sorting[, 21]),] -> temp.for_sorting

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

#bind all in one and leave only rows with the correct number of columns. So 3x3 + 3 + 3 + 1 = 16 or 3x5 + 3 + 5 + 1 = 24
do.call(rbind, do.call(rbind, reshaped.input.datalist)) -> temp.unbound
temp.unbound.filtered <- temp.unbound %>%
  dplyr::group_by(BlockID) %>%
    dplyr::mutate(nrows = n()) %>%
      dplyr::filter(nrows == 16 | nrows == 24)

#renames colnames
##tracemem says that by renaming columns we make it a completely new object. Which means copying over
##tracemem(unbound_filtered) -> before
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
do.call(rbind, ReshapeForDB(temp.reshaped.ZIP)) -> temp.rbound.ZIP
do.call(rbind, ReshapeForDB(temp.reshaped.Loewe)) -> temp.rbound.Loewe
do.call(rbind, ReshapeForDB(temp.reshaped.Bliss)) -> temp.rbound.Bliss
do.call(rbind, ReshapeForDB(temp.reshaped.HSA)) -> temp.rbound.HSA

temp.out <-  merge(merge(temp.rbound.ZIP, temp.rbound.Loewe, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID')), 
              merge(temp.rbound.HSA, temp.rbound.Bliss, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID')), 
              by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID'))

saveRDS(object = temp.out, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_All"))

#housekeeping
rm(list=ls(pattern='^temp'))
