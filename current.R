library("tidyverse")
library("synergyfinder")
setwd('/Users/zagidull/Documents/fimm_files/synergy_calc_august/almanac')


nov2017 <- read_csv('ComboDrugGrowth_Nov2017.csv', progress = F) # it is nice to read as tibble, this helps get rid of formatting errors
nov2017$TESTDATE <- as.Date(nov2017$TESTDATE,"%m/%d/%Y") #correcting DATA col class
sapply(nov2017, class)

#reading in
for_sorting <- nov2017
cellnames<-c('786-0','A498')

#two celllines
for_sorting_two <- for_sorting %>% # taking one cell line as an example
  filter(CELLNAME=='A498')
for_sorting_three <- for_sorting %>% # taking one cell line as an example/ NB used to be one with cell name 786-0.Did that not to change teh rest of the code
  filter(CELLNAME %in% cellnames)

#dropping columns that have one unique value. `dropped` var contains the names of these variables
before <- colnames(for_sorting_three)
for_sorting_three <- for_sorting_three[vapply(for_sorting_three, function(x) length(unique(x)) > 1, logical(1L))]
after <- colnames(for_sorting_three)
dropped <- setdiff(before, after)

for_sorting_one <- for_sorting_one[vapply(for_sorting_one, function(x) length(unique(x)) > 1, logical(1L))]
for_sorting_two <- for_sorting_two[vapply(for_sorting_two, function(x) length(unique(x)) > 1, logical(1L))]

#this is being run
library("tidyverse")
library("synergyfinder")
setwd('/Users/zagidull/Documents/fimm_files/synergy_calc_august/almanac')


nov2017 <- read_csv('ComboDrugGrowth_Nov2017.csv', progress = F) # it is nice to read as tibble, this helps get rid of formatting errors
nov2017$TESTDATE <- as.Date(nov2017$TESTDATE,"%m/%d/%Y") #correcting DATA col class
sapply(nov2017, class)

#reading in
for_sorting <- nov2017
for_sorting_one <- for_sorting %>% # taking one cell line as an example////nb this is the original input file
  filter(CELLNAME=='786-0')

# this filters out unique drug combos
uniqueDrug1Drub2Combos <- for_sorting_one %>%
  ungroup() %>%
  mutate(NSC1,NSC2) %>%
  select(NSC1, NSC2) %>%
  unique %>%
  na.omit()

#testing on first two
uniqueDrug1Drub2Combos[c(1:100),] -> test
uniqueDrug1Drub2Combos -> test


datalist_full <- list()
#x= for_sorting_one
y = test # exchange to uniqueDrug1Drug2Combos
for (i in 1:nrow(y)) {
  perrow <- list()
  input <- list()
  output1 <- list()
  param1 <- vector(mode="numeric", length=1)
  param2 <- vector(mode="numeric", length=1)
  perrow<- y[i,]
  
  param1 <- as.numeric(perrow[1])
  param2 <- as.numeric(perrow[2])
  
  input <- correctingConc(for_sorting_one, param1, param2)
  iden <- paste(i, "786-0", sep = ":")
  output1 <-  tryCatch(rowInsert(input,iden), error = function(e) NA)
  #    output1 <- rowInsert(input,i)
  datalist_full[[i]] <- output1
  
}

#datalist_full %>%
#  filter(length(datalist_full) != 16)

#this tests for a condition, when it is NSC2 tested as a single drug


# this exchanges concindex2 with concindex1, when predicate function returns true. Meaning 


# this adds the zero row and correct row columns to the input dataframe. NSC1 becomes column, NSC2 becomes row
rowInsert <- function(x, blockID){
  class_old <- list()
  class_old <- sapply(x, class)
  #  print(class_old)
  NSC2_un <- vector(mode="numeric", length=1)
  NSC2_un <- unique(x$NSC2[!is.na(x$NSC2)])
  predicate <- function(x) {
    if (x[["CONCINDEX2"]]== 0 & x[["NSC1"]] == NSC2_un) {
      return(T)
    }
    else {return(F)}
  }
  exchange <- function(x) {
    if (predicate(x)){
      placeholderIndex <- x[["CONCINDEX1"]]
      placeholderConc <- x[["CONC1"]]
      x[["CONCINDEX1"]] <- x[["CONCINDEX2"]]
      x[["CONCINDEX2"]] <- placeholderIndex
      x[["CONC1"]] <- x[["CONC2"]]
      x[["CONC2"]] <- placeholderConc
      return(x)
    } else { return(x)}
  } 
 
  x <- as.tibble(t(apply(x, 1, exchange)))
  x[] <- Map(`class<-`,x, class_old) #brilliant! https://stackoverflow.com/questions/27435873/changing-class-of-data-frame-columns-using-strings
  #  class_new <- sapply(x, class)
  #  print(class_new)
  
  
  
  
  NSC1 <- dplyr::setdiff(unique(x$NSC1[!is.na(x$NSC1)]),NSC2_un)
  #  NSC1 <- unique(x$NSC1[!is.na(x$NSC1)])
  x[(nrow(x)+1),] <- x[nrow(x),]
  x[nrow(x),]$NSC1 <- NSC1
  x[nrow(x),]$NSC2 <- NSC2_un
  x[nrow(x),]$CONCINDEX1 <- x[nrow(x),]$CONCINDEX2 <- x[nrow(x),]$CONC1 <- x[nrow(x),]$CONC2 <- 0
  x[nrow(x),]$mean_single <- 100
  x[nrow(x),]$SAMPLE2 <- unique(x$SAMPLE2[!is.na(x$SAMPLE2)])
  
  x$Row <- x$CONCINDEX2 + 1
  x$Col <- x$CONCINDEX1 + 1
  x$BlockID <- c(blockID) # sub with cell-line
  x$Response <- x$mean_single
  x$ConcRowUnit <- x$ConcColUnit <- c("M")
  x$ConcRow <- x$CONC2
  x$ConcCol <- x$CONC1
  x$ConcCol[is.na(x$ConcCol) ] <- 0
  x$DrugRow <- NSC2_un
  x$DrugCol <- NSC1
  x[is.na(x$ConcRow),]$ConcRow <- 0
  x <- x[,-c(1:16)]
  return(x)
}



#this is used to select correct concentrations from the list based on input data.frame, and NSC1 and NSC2
correctingConc <- function(x,y,z) {
  y <- as.numeric(y)
  z <- as.numeric(z)
  x <- filter(x, ((NSC1 == y & NSC2 == z) | (NSC1 == y & is.na(NSC2)) | (NSC1 == z & is.na(NSC2)) ))
  # the issue with this list was the forth filter param, it basically meant that it possible i would have two blocks of tests. Where z is main, or y is main
 # x <- filter(x, ((NSC1 == y & NSC2 == z) | (NSC1 == y & is.na(NSC2)) | (NSC1 == z & is.na(NSC2)) | (NSC1 == z & NSC2 ==y)))
  aCONC1 <- unique(filter(x, (!is.na(CONC2))  )$CONC1)
  aNSC1 <- y #unique(filter(x, (!is.na(CONC2))  )$NSC1)
  aCONC2 <- unique(filter(x, (!is.na(CONC2))  )$CONC2)
  aNSC2 <-z  # unique(filter(x, (!is.na(CONC2))  )$NSC2)
  x <- x %>%
    group_by(CONCINDEX1,CONC1, CONC2, NSC1, NSC2) %>%
    mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
    arrange()
  x <- filter(x, ( (!is.na(NSC1) & !is.na(NSC2)) | ((CONC1 %in% aCONC1)& NSC1 ==aNSC1) | ((CONC1 %in% aCONC2)&NSC1==aNSC2) )  )
  
  #here we operate under the assumption that the same drug should be tested by the same screening facility.So using prettyR::Mode()
  screenID <- prettyR::Mode(x$SCREENER)

  x <- x[!duplicated(x[c("NSC1","NSC2" ,"CONC1", "CONC2", "mean_single")]),]
  x <- filter(x,SCREENER ==  screenID)
  x <- x[-match(c("SCORE", "EXPECTEDGROWTH", "TZVALUE","CONTROLVALUE","TESTVALUE", "TESTDATE","STUDY"), names(x))]
  x <- x %>%
    filter(!(CONCINDEX1 < 0))
  
  return(x)
}

#looping over all the cellnames

datalist_full <- list()


cellnames<-unique(nov2017$CELLNAME)

for (name in cellnames) {
  #this creates a place holder for one cellline 
  for_sorting_many <- list()
  for_sorting_many <- for_sorting %>% # taking one cell line as an example////nb this is the original input file
    filter(CELLNAME==name)
  
  #this makes a list of drugs used in given cellline
  uniqueDrug1Drub2Combos <- list()
  uniqueDrug1Drub2Combos <- for_sorting_many %>%
    mutate(NSC1,NSC2) %>%
    select(NSC1, NSC2) %>%
    unique %>%
    na.omit()
  
  y = uniqueDrug1Drub2Combos # exchange to uniqueDrug1Drug2Combos
  for (i in 1:nrow(y)) {
    perrow <- list()
    input <- list()
    output1 <- list()
    param1 <- vector(mode="numeric", length=1)
    param2 <- vector(mode="numeric", length=1)
    perrow<- y[i,]
    
    param1 <- as.numeric(perrow[1])
    param2 <- as.numeric(perrow[2])
    
    input <- correctingConc(for_sorting_one, param1, param2)
    iden <- paste(i, name, sep = ":")
    output1 <-  tryCatch(rowInsert(input,iden), error = function(e) NA)
    #    output1 <- rowInsert(input,i)
    datalist_full[[name]][[i]] <- output1
    
  }
  
  
  }



# for when concentrations tested of single drug are not present in combos. 
test2 <- test %>%
  filter(!(CONCINDEX1 < 0))

# for timebeing separation to leave only two drugs
test2 <- test2 %>%
  filter(!(NSC1 == 256942))



input <- filter(for_sorting_one, ((NSC1 == 256942 & NSC2 == 707389) | (NSC1 == 707389 & is.na(NSC2)) | (NSC1 == 256942 & is.na(NSC2)) | (NSC1 == 707389 & NSC2 ==256942)))
aCONC1 <- unique(filter(input, (!is.na(CONC2))  )$CONC1)
aNSC1 <- unique(filter(input, (!is.na(CONC2))  )$NSC1)
aCONC2 <- unique(filter(input, (!is.na(CONC2))  )$CONC2)
aNSC2 <- unique(filter(input, (!is.na(CONC2))  )$NSC2)
input <- filter(input, ( (!is.na(NSC1) & !is.na(NSC2)) | ((CONC1 %in% aCONC1)& NSC1 ==aNSC1) | ((CONC1 %in% aCONC2)&NSC1==aNSC2) )  )
input <- input %>%
  group_by(CONCINDEX1,CONC1, CONC2, NSC1, NSC2) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()

input <- input[!duplicated(input[c("NSC1","NSC2" ,"CONC1", "CONC2", "mean_single")]),]

input <- input %>%
  filter(!(CONCINDEX1 < 0))


# for filtering based on correct CONC and NSC values
filter(input, ( (!is.na(NSC1) & !is.na(NSC2)) | ((CONC1 %in% aCONC1)& NSC1 ==aNSC1) | ((CONC1 %in% aCONC2)&NSC1==aNSC2) )  )

aCONC1 <- unique(filter(input, (!is.na(CONC2))  )$CONC1)
aNSC1 <- unique(filter(input, (!is.na(CONC2))  )$NSC1)
aCONC2 <- unique(filter(input, (!is.na(CONC2))  )$CONC2)
aNSC2 <- unique(filter(input, (!is.na(CONC2))  )$NSC2)

#dropping drugs that have only been tested as a single. But wait, do i really need to drop it? I do because 
singles <- for_sorting_one %>% # this is for al lthe cases when only a single drug was tested
  filter(CONCINDEX2 == 0)
combos <- for_sorting_one %>% # this is for al lthe cases when two drugs were tested
  filter(CONCINDEX2 != 0)
NSC_to_drop <- setdiff(singles$NSC1, combos$NSC1)
singles <- singles[-which(singles$NSC1 == NSC_to_drop),]
for_syn_calc <- rbind(singles, combos)



#adding groupID index. For usage in pandas
for_sorting_three$groupID <- for_sorting_three %>% 
  group_by(CELLNAME) %>%
  group_indices(.)
for_sorting_three <- for_sorting_three %>%
  group_by(groupID)



