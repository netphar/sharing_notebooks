# MIT license
# this file uses NCI ALMANAC data as input: 
# https://wiki.nci.nih.gov/download/attachments/338237347/ComboDrugGrowth_Nov2017.zip?version=1&modificationDate=1510057275000&api=v2
# it reshapes and modifies to be used fpr synergyfinder::ReshapeData()
# authors: Tang, J., Zagidullin, B.
# version 0.3

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

#input file
file <- c('ComboDrugGrowth_Nov2017.csv')
#file <- c('new.csv')

#reading in the data and correcting col classes
input.data <- read_csv(file, progress = F) # it is nice to read as tibble, this helps get rid of formatting errors
#input.data.full <- read_csv(file1, progress = F) # it is nice to read as tibble, this helps get rid of formatting errors

input.data$TESTDATE <- as.Date(input.data$TESTDATE,"%m/%d/%Y") #correcting DATA col class
#input.data.full$TESTDATE <- as.Date(input.data$TESTDATE,"%m/%d/%Y")
#sapply(temp.for_sorting, class) this checks the class of each of the columns
input.data -> temp.for_sorting

#creating placeholde for results  
reshaped.input.datalist <- list()

#cellnames to be searched for
temp.cellnames_old <- unique(temp.for_sorting$CELLNAME)
temp.cellnames <- temp.cellnames_old[-61] #todo check why the last cell line has to be removed
#temp.cellnames <- temp.cellnames_old

#remove rows where PERCENTGROWTHNOTZ is na. 
# https://stackoverflow.com/questions/4862178/remove-rows-with-all-or-some-nas-missing-values-in-data-frame
temp.for_sorting[complete.cases(temp.for_sorting[, 21]),] -> temp.for_sorting

#cellnames <- c("HCT-116")

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

#bind all in one and leave only rows with the correct number of columns. So 3x3 + 3 + 3 + 1 or 3x5 + 3 + 5 + 1
do.call(rbind, do.call(rbind, reshaped.input.datalist)) -> temp.unbound
temp.unbound.filtered <- temp.unbound %>%
  dplyr::group_by(BlockID) %>%
    dplyr::mutate(nrows = n()) %>%
      dplyr::filter(nrows == 16 | nrows == 24)

#rename colnames
##tracemem(unbound_filtered) -> before
names(temp.unbound.filtered) <- c('row', 'col', 'block_id', 'response', 'conc_c_unit', 'conc_r_unit', 'conc_r', 'conc_c', 'drug_row', 'drug_col', 'nrows')
##tracemem(unbound_filtered) -> after

#reshaping
reshaped <- ReshapeData(temp.unbound.filtered, data.type = 'viability')

#change concetrations from M to uM
reshaped$dose.response.mats <- lapply(reshaped$dose.response.mats, function(x) {
  colnames(x) <- as.numeric(colnames(x))*(10^6)
  rownames(x) <- as.numeric(rownames(x))*(10^6)
  return(x)
})
reshaped$drug.pairs$concRUnit <- reshaped$drug.pairs$concCUnit <- 'uM'

CalculateSynergy(reshaped, method = 'ZIP', correction = T, Emin = 0, Emax = NA) -> reshaped.ZIP
CalculateSynergy(reshaped, method = 'Loewe', correction = T, Emin = 0, Emax = NA) -> reshaped.Loewe
CalculateSynergy(reshaped, method = 'Bliss', correction = T, Emin = 0) -> reshaped.Bliss
CalculateSynergy(reshaped, method = 'HSA', correction = T, Emin = 0, Emax = NA) -> reshaped.HSA

# datalist
calculated.synergy.datalist = list()
temp.mylist <- list()
temp.mylist.Bliss <- list()
temp.mylist.HSA <- list()
temp.mylist.Loewe <- list()

temp.first.two <- list()
temp.first.two.Bliss <- list()
temp.first.two.HSA <- list()
temp.first.two.Loewe <- list()

#populate for ZIP
temp.first.two$dose.response.mats <- reshaped.ZIP$dose.response.mats
temp.first.two$drug.pairs <- reshaped.ZIP$drug.pairs
temp.first.two$scores <- reshaped.ZIP$scores
temp.first.two$method <- reshaped.ZIP$method

#populate for Bliss
temp.first.two.Bliss$dose.response.mats <- reshaped.Bliss$dose.response.mats
temp.first.two.Bliss$drug.pairs <- reshaped.Bliss$drug.pairs
temp.first.two.Bliss$scores <- reshaped.Bliss$scores
temp.first.two.Bliss$method <- reshaped.Bliss$method

#populate for HSA
temp.first.two.HSA$dose.response.mats <- reshaped.HSA$dose.response.mats
temp.first.two.HSA$drug.pairs <- reshaped.HSA$drug.pairs
temp.first.two.HSA$scores <- reshaped.HSA$scores
temp.first.two.HSA$method <- reshaped.HSA$method

#populate for Loewe
temp.first.two.Loewe$dose.response.mats <- reshaped.Loewe$dose.response.mats
temp.first.two.Loewe$drug.pairs <- reshaped.Loewe$drug.pairs
temp.first.two.Loewe$scores <- reshaped.Loewe$scores
temp.first.two.Loewe$method <- reshaped.Loewe$method

for (i in 1:nrow(temp.first.two$drug.pairs)) {
  temp.mylist[i] <- temp.first.two$drug.pairs[i,]$blockIDs
}

for (i in 1:nrow(temp.first.two.Bliss$drug.pairs)) {
  temp.mylist.Bliss[i] <- temp.first.two.Bliss$drug.pairs[i,]$blockIDs
}
for (i in 1:nrow(temp.first.two.HSA$drug.pairs)) {
  temp.mylist.HSA[i] <- temp.first.two.HSA$drug.pairs[i,]$blockIDs
}

for (i in 1:nrow(temp.first.two.Loewe$drug.pairs)) {
  temp.mylist.Loewe[i] <- temp.first.two.Loewe$drug.pairs[i,]$blockIDs
}

temp.pb <- txtProgressBar(min = 0, max = nrow(temp.first.two$drug.pairs), style = 3)

for (i in 1:nrow(temp.first.two$drug.pairs))
{
  setTxtProgressBar(temp.pb, i)
  temp.a <- list()
  temp.b <- list()
  temp.c <- list()
  temp.d <- list()
  
  
  temp.a$dose.response.mats <- temp.first.two$dose.response.mats[i]
  temp.a$scores <- temp.first.two$scores[i]
  temp.a$method <- temp.first.two$method
  temp.a$drug.pairs <- temp.first.two$drug.pairs[i,]
  
  if (is.na(temp.a$scores)) {
    cols <- colnames(temp.a$dose.response.mats[[which(is.na(temp.a$scores))]])
    rows <- rownames(temp.a$dose.response.mats[[which(is.na(temp.a$scores))]])
    numrows <- dim(temp.a$dose.response.mats[[which(is.na(temp.a$scores))]])[1]
    numcols <- dim(temp.a$dose.response.mats[[which(is.na(temp.a$scores))]])[2]
    temp.a$scores[[which(is.na(temp.a$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  temp.b$dose.response.mats <- temp.first.two.Bliss$dose.response.mats[i]
  temp.b$scores <- temp.first.two.Bliss$scores[i]
  temp.b$method <- temp.first.two.Bliss$method
  temp.b$drug.pairs <- temp.first.two.Bliss$drug.pairs[i,]
  
  if (is.na(temp.b$scores)) {
    cols <- colnames(temp.b$dose.response.mats[[which(is.na(temp.b$scores))]])
    rows <- rownames(temp.b$dose.response.mats[[which(is.na(temp.b$scores))]])
    numrows <- dim(temp.b$dose.response.mats[[which(is.na(temp.b$scores))]])[1]
    numcols <- dim(temp.b$dose.response.mats[[which(is.na(temp.b$scores))]])[2] 
    temp.b$scores[[which(is.na(temp.b$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  temp.c$dose.response.mats <- temp.first.two.Loewe$dose.response.mats[i]
  temp.c$scores <- temp.first.two.Loewe$scores[i]
  temp.c$method <- temp.first.two.Loewe$method
  temp.c$drug.pairs <- temp.first.two.Loewe$drug.pairs[i,]
  
  if (is.na(temp.c$scores)) {
    
    cols <- colnames(temp.c$dose.response.mats[[which(is.na(temp.c$scores))]])
    rows <- rownames(temp.c$dose.response.mats[[which(is.na(temp.c$scores))]])
    numrows <- dim(temp.c$dose.response.mats[[which(is.na(temp.c$scores))]])[1]
    numcols <- dim(temp.c$dose.response.mats[[which(is.na(temp.c$scores))]])[2]
    temp.c$scores[[which(is.na(temp.c$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  temp.d$dose.response.mats <- temp.first.two.HSA$dose.response.mats[i]
  temp.d$scores <- temp.first.two.HSA$scores[i]
  temp.d$method <- temp.first.two.HSA$method
  temp.d$drug.pairs <- temp.first.two.HSA$drug.pairs[i,]
  
  if (is.na(temp.d$scores)) {
    cols <- colnames(temp.d$dose.response.mats[[which(is.na(temp.d$scores))]])
    rows <- rownames(temp.d$dose.response.mats[[which(is.na(temp.d$scores))]])
    numrows <- dim(temp.d$dose.response.mats[[which(is.na(temp.d$scores))]])[1]
    numcols <- dim(temp.d$dose.response.mats[[which(is.na(temp.d$scores))]])[2]
    temp.d$scores[[which(is.na(temp.d$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  temp.x.ZIP <- melt(temp.a$dose.response.mats)
  temp.y.ZIP <- melt(temp.a$scores)
  temp.xy.ZIP <- merge(temp.x.ZIP,temp.y.ZIP, by = c('Var1', 'Var2', 'L1'))
  
  temp.x.Bliss <- melt(temp.b$dose.response.mats)
  temp.y.Bliss <- melt(temp.b$scores)
  temp.xy.Bliss <- merge(temp.x.Bliss,temp.y.Bliss, by = c('Var1', 'Var2', 'L1'))
  
  temp.x.Loewe <- melt(temp.c$dose.response.mats)
  temp.y.Loewe <- melt(temp.c$scores)
  temp.xy.Loewe <- merge(temp.x.Loewe,temp.y.Loewe, by = c('Var1', 'Var2', 'L1'))
  
  temp.x.HSA <- melt(temp.d$dose.response.mats)
  temp.y.HSA <- melt(temp.d$scores)
  temp.xy.HSA <- merge(temp.x.HSA,temp.y.HSA, by = c('Var1', 'Var2', 'L1'))
  
  colnames(temp.xy.ZIP)[which(names(temp.xy.ZIP) == "Var1")] <- "ConcR"
  colnames(temp.xy.ZIP)[which(names(temp.xy.ZIP) == "Var2")] <- "ConcC"
  colnames(temp.xy.ZIP)[which(names(temp.xy.ZIP) == "L1")] <- "blockIDs"
  temp.xy.ZIP$blockIDs <- temp.mylist[[i]]
  colnames(temp.xy.ZIP)[which(names(temp.xy.ZIP) == "value.x")] <- "Response_inhibition"
  temp.synergy.type.ZIP <- paste('Synergy', temp.first.two$method, sep = '_')
  colnames(temp.xy.ZIP)[which(names(temp.xy.ZIP) == "value.y")] <- temp.synergy.type.ZIP
  
  colnames(temp.xy.Bliss)[which(names(temp.xy.Bliss) == "Var1")] <- "ConcR"
  colnames(temp.xy.Bliss)[which(names(temp.xy.Bliss) == "Var2")] <- "ConcC"
  colnames(temp.xy.Bliss)[which(names(temp.xy.Bliss) == "L1")] <- "blockIDs"
  temp.xy.Bliss$blockIDs <- temp.mylist.Bliss[[i]]
  colnames(temp.xy.Bliss)[which(names(temp.xy.Bliss) == "value.x")] <- "Response_inhibition"
  temp.synergy.type.Bliss <- paste('Synergy', temp.first.two.Bliss$method, sep = '_')
  colnames(temp.xy.Bliss)[which(names(temp.xy.Bliss) == "value.y")] <- temp.synergy.type.Bliss
  
  colnames(temp.xy.Loewe)[which(names(temp.xy.Loewe) == "Var1")] <- "ConcR"
  colnames(temp.xy.Loewe)[which(names(temp.xy.Loewe) == "Var2")] <- "ConcC"
  colnames(temp.xy.Loewe)[which(names(temp.xy.Loewe) == "L1")] <- "blockIDs"
  temp.xy.Loewe$blockIDs <- temp.mylist.Loewe[[i]]
  colnames(temp.xy.Loewe)[which(names(temp.xy.Loewe) == "value.x")] <- "Response_inhibition"
  temp.synergy.type.Loewe <- paste('Synergy', temp.first.two.Loewe$method, sep = '_')
  colnames(temp.xy.Loewe)[which(names(temp.xy.Loewe) == "value.y")] <- temp.synergy.type.Loewe
  
  colnames(temp.xy.HSA)[which(names(temp.xy.HSA) == "Var1")] <- "ConcR"
  colnames(temp.xy.HSA)[which(names(temp.xy.HSA) == "Var2")] <- "ConcC"
  colnames(temp.xy.HSA)[which(names(temp.xy.HSA) == "L1")] <- "blockIDs"
  temp.xy.HSA$blockIDs <- temp.mylist.HSA[[i]]
  colnames(temp.xy.HSA)[which(names(temp.xy.HSA) == "value.x")] <- "Response_inhibition"
  temp.synergy.type.HSA <- paste('Synergy', temp.first.two.HSA$method, sep = '_')
  colnames(temp.xy.HSA)[which(names(temp.xy.HSA) == "value.y")] <- temp.synergy.type.HSA
  
  
  temp.first.two.no.celllines <- merge(temp.xy.ZIP, temp.a$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  temp.first.two.no.celllines.Bliss <- merge(temp.xy.Bliss, temp.b$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  temp.first.two.no.celllines.Loewe <- merge(temp.xy.Loewe, temp.c$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  temp.first.two.no.celllines.HSA <- merge(temp.xy.HSA, temp.d$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  
  temp.first.two.no.celllines.ZIP.Bliss <- merge(temp.first.two.no.celllines, temp.first.two.no.celllines.Bliss, 
                                            by = c('blockIDs','ConcR','ConcC','Response_inhibition',
                                                   'drug.row', 'drug.col','concRUnit','concCUnit'))
  temp.first.two.no.celllines.ZIP.Bliss.HSA <- merge(temp.first.two.no.celllines.ZIP.Bliss, temp.first.two.no.celllines.HSA, 
                                                by = c('blockIDs','ConcR','ConcC','Response_inhibition',
                                                       'drug.row', 'drug.col','concRUnit','concCUnit'))
  temp.first.two.no.celllines.ZIP.Bliss.HSA.Loewe <- merge(temp.first.two.no.celllines.ZIP.Bliss.HSA, temp.first.two.no.celllines.Loewe, 
                                                      by = c('blockIDs','ConcR','ConcC','Response_inhibition',
                                                             'drug.row', 'drug.col','concRUnit','concCUnit'))
  
calculated.synergy.datalist[[i]] <- temp.first.two.no.celllines.ZIP.Bliss.HSA.Loewe
  
}
close(temp.pb)
calculated.synergy.rbound.datalist = do.call(rbind, calculated.synergy.datalist)
#saveRDS(object = big_data, file = '3008_big_data_with_Loewe') #result in table format
#saveRDS(object = datalist, file = '3008_datalist_with_Loewe') # result as nested list


#housekeeping
rm(list=ls(pattern='^temp'))
