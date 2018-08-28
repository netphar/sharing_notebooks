#housekeeping
library("tidyverse")
library('drc')
library('reshape2') # is it necessary?
#library("synergyfinder") # i do not want to use save
setwd('/Users/zagidull/Documents/fimm_files/synergy_calc_august/almanac')

#reading in the data and correcting col classes
nov2017 <- read_csv('ComboDrugGrowth_Nov2017.csv', progress = F) # it is nice to read as tibble, this helps get rid of formatting errors
nov2017$TESTDATE <- as.Date(nov2017$TESTDATE,"%m/%d/%Y") #correcting DATA col class
sapply(nov2017, class)
nov2017 -> for_sorting

# this adds the zero row and correct row columns to the input dataframe. NSC1 becomes column, NSC2 becomes row
rowInsert <- function(x, blockID){
  class_old <- list()
  class_old <- sapply(x, class)
  #  print(class_old)
  NSC2_un <- vector(mode="numeric", length=1)
  NSC2_un <- unique(x$NSC2[!is.na(x$NSC2)])

  predicate <- function(x) {
    if (as.numeric(x[["CONCINDEX2"]]== 0) & as.numeric(x[["NSC1"]]) == NSC2_un) {
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
#  x <- x[,-c(1:16)] # to-do change that
  x <- x[-match(c("COMBODRUGSEQ","SCREENER","STUDY", "TESTDATE", "PLATE", "PANELNBR", "CELLNBR", "PREFIX1", "NSC1", "SAMPLE1", "CONCINDEX1", "CONC1", "CONCUNIT1", "PREFIX2", "NSC2", "SAMPLE2", "CONCINDEX2", "CONC2", "CONCUNIT2", "PERCENTGROWTH", "PERCENTGROWTHNOTZ","TESTVALUE", "CONTROLVALUE", "TZVALUE", "EXPECTEDGROWTH", "SCORE", "VALID", "PANEL", "CELLNAME", "mean_single"), names(x))]
  return(x)
}

#this is used to select correct concentration and drug combinations
correctingConc <- function(x,y,z) {
  
  y <- as.numeric(y)
  z <- as.numeric(z)

  x <- filter(x, ((NSC1 == y & NSC2 == z) | (NSC1 == y & is.na(NSC2)) | (NSC1 == z & is.na(NSC2)) ))
  
  x_comb <- list()
  x_comb <- filter(x, ( (NSC1 == y & NSC2 == z) ) )
  x_comb <- x_comb %>%
    group_by(NSC1, NSC2, CONC1, CONC2) %>%
    mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
    arrange()
  
  
  aCONC11 <- unique(filter(x_comb, (!is.na(CONC2))  )$CONC1)
  aNSC11 <- y
  aCONC21 <- unique(filter(x_comb, (!is.na(CONC2))  )$CONC2)
  aNSC21 <-z
  
  x_one <- filter(x, (CONC1 %in% aCONC11)& NSC1 ==aNSC11 & is.na(CONC2))
  x_one <- x_one %>%
    group_by(CONC1) %>%
    mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
    arrange()
  df <- data.frame()
  df <- as.data.frame(cbind(order(aCONC11), aCONC11))
  x_one <- inner_join(x_one, df, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC11"))
  x_one <- x_one[!duplicated(x_one[c("CONC1")]),]
  
  x_two <- filter(x, (CONC1 %in% aCONC21)& NSC1 ==aNSC21 & is.na(CONC2))
  x_two <- x_two %>%
    group_by(CONC1) %>%
    mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
    arrange()
  df2 <- data.frame()
  df2 <- as.data.frame(cbind(order(aCONC21), aCONC21))
  x_two <- inner_join(x_two, df2, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC21"))
  x_two <- x_two[!duplicated(x_two[c("CONC1")]),]
  
  x_comb <- rbind(x_comb, x_one)
  x_comb <- rbind(x_comb, x_two)
  return(x_comb)
}

#creating placeholde for results  
datalist_full <- list()

#cellnames to be searched for
cellnames_old <- unique(nov2017$CELLNAME)
cellnames <- cellnames_old[-61] #todo check why the last cell line has to be removed

#cellnames <- c("HCT-116")

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
    
    input <- correctingConc(for_sorting_many, param1, param2)
    iden <- paste(i, name, sep = ":")
    output1 <-  tryCatch(rowInsert(input,iden), error = function(e) NA)
    #    output1 <- rowInsert(input,i)
    datalist_full[[name]][[i]] <- output1
    if (i%%1000==0) cat('left', nrow(y) - i,'\n',sep='' )
    if (i == nrow(y)) cat('cell group processed = ', name, ' \n', sep='')
    
  }
  
}

#not sure what s the better way or rather if there is any difference between thes two. Choose either. resuts and runtime are identical
# one
do.call(rbind, datalist_full) -> one
do.call(rbind, one) -> two
two_cleaned <- two %>% group_by(BlockID) %>% mutate(nrows = n()) %>% filter(nrows == 16 | nrows == 24)
# two
do.call(rbind, do.call(rbind, datalist_full)) -> unbound
unbound_filtered <- unbound %>%
  group_by(BlockID) %>%
    mutate(nrows = n()) %>%
      filter(nrows == 16 | nrows == 24)

#reshaping
reshaped <- ReshapeData(unbound_filtered, data.type = 'viability')


#this calculates synergy for the full valid 3x3 and 3x5 matrixes of NCI ALMANAC. 307737 dose response mats
# ZIP, HSA and Bliss were calculated on macosx with the following sessioninfo()
#########################################################################################################################################################
#R version 3.4.4 (2018-03-15)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS High Sierra 10.13.6
#
#Matrix products: default
#BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
#
#locale:
#  [1] C
#
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#  [1] synergyfinder_1.4.2 nleqslv_3.3.2       rlist_0.4.6.1       drc_3.0-1           MASS_7.3-50         reshape2_1.4.3      svMisc_1.1.0        prettyR_2.2         bindrcpp_0.2.2     
#[10] forcats_0.3.0       stringr_1.3.1       dplyr_0.7.5         purrr_0.2.5         readr_1.1.1         tidyr_0.8.1         tibble_1.4.2        ggplot2_2.2.1       tidyverse_1.2.1    
#
#loaded via a namespace (and not attached):
#  [1] httr_1.3.1            maps_3.3.0            jsonlite_1.5          splines_3.4.4         carData_3.0-1         modelr_0.1.2          gtools_3.5.0          assertthat_0.2.0     
#[9] cellranger_1.1.0      yaml_2.1.19           pillar_1.2.3          lattice_0.20-35       glue_1.2.0            rvest_0.3.2           colorspace_1.3-2      sandwich_2.4-0       
#[17] Matrix_1.2-14         plyr_1.8.4            psych_1.8.4           pkgconfig_2.0.1       broom_0.4.4           haven_1.1.1           mvtnorm_1.0-8         scales_0.5.0         
#[25] gdata_2.18.0          openxlsx_4.1.0        rio_0.5.10            car_3.0-0             TH.data_1.0-8         lazyeval_0.2.1        cli_1.0.0             mnormt_1.5-5         
#[33] survival_2.42-3       magrittr_1.5          crayon_1.3.4          readxl_1.1.0          nlme_3.1-137          gplots_3.0.1          xml2_1.2.0            foreign_0.8-70       
#[41] tools_3.4.4           SpatialExtremes_2.0-6 data.table_1.11.4     hms_0.4.2             multcomp_1.4-8        gridBase_0.4-7        munsell_0.5.0         plotrix_3.7-2        
#[49] zip_1.0.0             compiler_3.4.4        caTools_1.17.1        rlang_0.2.1           grid_3.4.4            rstudioapi_0.7        bitops_1.0-6          gtable_0.2.0         
#[57] codetools_0.2-15      abind_1.4-5           curl_3.2              R6_2.2.2              zoo_1.8-2             lubridate_1.7.4       utf8_1.1.4            bindr_0.1.1          
#[65] KernSmooth_2.23-15    stringi_1.2.3         parallel_3.4.4        Rcpp_0.12.17          tidyselect_0.2.4    
#########################################################################################################################################################


CalculateSynergy(reshaped, method = "HSA") -> reshaped.HSA
CalculateSynergy(reshaped, method = "Bliss") -> reshaped.Bliss
CalculateSynergy(reshaped, method = "ZIP", correction = T) -> reshaped.ZIP
#save using saveRDS()


# for synergy calc
#library("drc")



#it tries to use only LL.4 model for fitting the data. 
# but d1.fun and d2.fun depend one which model was used to fit: LL.4 vs L.4
# so we end up having 4 scenarios -> both row and col LL4, both row and col L4, and then row LL4 col L4; row L4 and col LL4
# drug.col.par/drug.row.par 1,2,3,4 -> b,c,d,e

#original
Loewe_modded_original <- function (response.mat, correction = T, Emin = 0, Emax = 100, 
                          nan.handle = c("LL4", "L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD(response.mat, Emin = Emin, 
                                         Emax = Emax, nan.handle)$corrected.mat
  }
  single.fit <- FittingSingleDrug(response.mat)
  drug.col.model <- single.fit$drug.col.model
  drug.col.par <- coef(drug.col.model)
  d1.fun <- function(conc, drug.col.par) {
    (drug.col.par[3] + drug.col.par[2] * (conc/drug.col.par[4])^drug.col.par[1])/(1 + 
                                                                                    (conc/drug.col.par[4])^drug.col.par[1])
  }
  drug.row.model <- single.fit$drug.row.model
  drug.row.par <- coef(drug.row.model)
  d2.fun <- function(conc, drug.row.par) {
    (drug.row.par[3] + drug.row.par[2] * (conc/drug.row.par[4])^drug.row.par[1])/(1 + 
                                                                                    (conc/drug.row.par[4])^drug.row.par[1])
  }
  row.conc <- as.numeric(rownames(response.mat))[-1]
  col.conc <- as.numeric(colnames(response.mat))[-1]
  loewe.mat <- response.mat
  for (i in 1:length(col.conc)) {
    for (j in 1:length(row.conc)) {
      x1 <- col.conc[i]
      x2 <- row.conc[j]
      eq <- function(x) {
        x1/(drug.col.par[4] * (((x - drug.col.par[3])/(drug.col.par[2] - 
                                                         x))^(1/drug.col.par[1]))) + x2/(drug.row.par[4] * 
                                                                                           (((x - drug.row.par[3])/(drug.row.par[2] - 
                                                                                                                      x))^(1/drug.row.par[1]))) - 1
      }
      slv <- nleqslv(max(drug.col.par[2] + 1, drug.row.par[2] + 
                           1), eq, method = "Newton")
      if (slv$termcd == 1) {
        loewe.mat[j + 1, i + 1] <- slv$x
      }
      else {
        y.loewe1 <- d1.fun(x1 + x2, drug.col.par)
        y.loewe2 <- d2.fun(x1 + x2, drug.row.par)
        loewe.mat[j + 1, i + 1] <- ifelse(y.loewe1 > 
                                            y.loewe2, y.loewe1, y.loewe2)
      }
    }
  }
  return(response.mat - loewe.mat)
}

#modified by jing 27/08
Loewe_modded <- function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
                                nan.handle = c("L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD2(response.mat, Emin = Emin, 
                                          Emax = Emax, nan.handle)$corrected.mat
  }
  single.fit <- FittingSingleDrug2(response.mat,fixed = c(NA, NA, NA, NA), nan.handle)
  drug.col.model <- single.fit$drug.col.model
  drug.col.par <- coef(drug.col.model)
  d1.fun <- function(conc, drug.col.model) {
    drug.col.par <- coef(drug.col.model)
    # LL.4
    if(length(grep("LL.4", drug.col.model$call$fct))> 0 )
      (drug.col.par[3] + drug.col.par[2] * (conc/drug.col.par[4])^drug.col.par[1])/(1 + 
                                                                                      (conc/drug.col.par[4])^drug.col.par[1])
    else # L.4
      (drug.col.par[2] + (drug.col.par[3] - drug.col.par[2])/(1+exp(drug.col.par[1]*(conc-drug.col.par[4]))))
  }
  
  drug.row.model <- single.fit$drug.row.model
  drug.row.par <- coef(drug.row.model)
  d2.fun <- function(conc, drug.row.model) {
    drug.row.par <- coef(drug.row.model)
    if(length(grep("LL.4", drug.row.model$call$fct))> 0 ) # LL.4
      (drug.row.par[3] + drug.row.par[2] * (conc/drug.row.par[4])^drug.row.par[1])/(1 + 
                                                                                      (conc/drug.row.par[4])^drug.row.par[1])
    else #L.4
      (drug.row.par[2] + (drug.row.par[3] - drug.row.par[2])/(1+exp(drug.row.par[1]*(conc-drug.row.par[4]))))
  }
  row.conc <- as.numeric(rownames(response.mat))[-1]
  col.conc <- as.numeric(colnames(response.mat))[-1]
  loewe.mat <- response.mat
  
  # Four functions to calculate loewe 
  eq.LL4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {# Eq.8 in the ZIP paper
    x1/(drug.col.par[4] * (((x - drug.col.par[3])/(drug.col.par[2] - x))^(1/drug.col.par[1]))) + 
      x2/(drug.row.par[4] * (((x - drug.row.par[3])/(drug.row.par[2] - x))^(1/drug.row.par[1]))) - 1
  }
  
  eq.L4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
    x1/(drug.col.par[4] + log((drug.col.par[3]-x)/(x-drug.col.par[2]))/drug.col.par[1]) +
      x2/(drug.row.par[4] + log((drug.row.par[3]-x)/(x-drug.row.par[2]))/drug.row.par[1]) -1
  }
  
  eq.LL4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
    x1/(drug.col.par[4] * (((x - drug.col.par[3])/(drug.col.par[2] - x))^(1/drug.col.par[1]))) +
      x2/(drug.row.par[4] + log((drug.row.par[3]-x)/(x-drug.row.par[2]))/drug.row.par[1]) -1
  }
  
  eq.L4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
    x1/(drug.col.par[4] + log((drug.col.par[3]-x)/(x-drug.col.par[2]))/drug.col.par[1]) +
      x2/(drug.row.par[4] * (((x - drug.row.par[3])/(drug.row.par[2] - x))^(1/drug.row.par[1]))) - 1
  }
  
  cond1 = length(grep("LL.4", drug.col.model$call$fct))> 0 
  cond2 = length(grep("LL.4", drug.row.model$call$fct))> 0 
  if( cond1 == T & cond2 == T) eq = eq.LL4.LL4
  if( cond1 == T & cond2 == F) eq = eq.LL4.L4
  if( cond1 == F & cond2 == T) eq = eq.L4.LL4
  if( cond1 == F & cond2 == F) eq = eq.L4.L4
  
  for (i in 1:length(col.conc)) {
    for (j in 1:length(row.conc)) {
      x1 <- col.conc[i]
      x2 <- row.conc[j]
      
      options(warn = -1)
      slv = tryCatch(
        {
          slv <- nleqslv(max(drug.col.par[2] + 1, drug.row.par[2] + 1), eq, method = "Newton", x1=x1, x2=x2, drug.col.par = drug.col.par, drug.row.par = drug.row.par)
          
        }, error = function(cond){
          slv = list(termcd = 0, x = 0)
        }
      )
      
      
      options(warn = 0)
      
      if (slv$termcd == 1) {
        loewe.mat[j + 1, i + 1] <- slv$x
      }
      else {
        y.loewe1 <- d1.fun(x1 + x2, drug.col.model)
        y.loewe2 <- d2.fun(x1 + x2, drug.row.model)
        loewe.mat[j + 1, i + 1] <- ifelse(100 > max(y.loewe1, y.loewe2), max(y.loewe1, y.loewe2), 100)
      }
    }
  }
  
  return(response.mat - loewe.mat)
}
Loewe <- function (response.mat,correction = T, Emin = 0, Emax = NA,nan.handle = c("L4")) {
  # scores <- list()
  #  method <- "Loewe"
  out <- tryCatch(Loewe_modded(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
                               nan.handle = nan.handle), error = function(e) NA)
  #  {
  #  x$method <- method
  #  x$scores <- x$dose.response.mats
  #  x$scores[!is.na(x$scores)] <- NA
  #  x$scores <- replace(x$scores, !is.na(x$scores), NA)
  #  purrr::map(!is.na(x$scores), NA)
  #  x$scores <- lapply(x$scores, function(z) ifelse(!is.na(z), NA, z))
  #  return(x)
  #  })
  return(out)
}

ZIP_modded <- function (response.mat, correction = T, Emin = 0, Emax = 100, 
                        nan.handle = c("LL4", "L4")) 
{
  if (correction) {
    nan.handle <- match.arg(nan.handle)
    response.mat <- BaselineCorrectionSD(response.mat, NA, 
                                         NA, nan.handle)$corrected.mat
  }
  single.fitted <- FittingSingleDrug(response.mat, fixed = c(NA, 
                                                             Emin, Emax, NA), nan.handle)
  drug.col.response <- single.fitted$drug.col.fitted
  drug.row.response <- single.fitted$drug.row.fitted
  updated.single.mat <- mat.or.vec(nrow(response.mat), ncol(response.mat))
  colnames(updated.single.mat) <- colnames(response.mat)
  rownames(updated.single.mat) <- rownames(response.mat)
  #oringinal
  #updated.single.mat[1, c(2:ncol(response.mat))] <- drug.col.response
  #updated.single.mat[c(2:nrow(response.mat)), 1] <- drug.row.response
  #modified by Jing 27/08
  updated.single.mat[1, c(2:ncol(response.mat))] <- drug.col.response[-1]
  updated.single.mat[c(2:nrow(response.mat)), 1] <- drug.row.response[-1]
  
  updated.col.mat <- updated.single.mat
  
  for (i in 2:ncol(response.mat)) {
    tmp <- as.data.frame(mat.or.vec(nrow(response.mat) - 
                                      1, 0))
    tmp$dose <- as.numeric(rownames(response.mat)[-1])
    tmp$inhibition <- response.mat[c(2:nrow(response.mat)), 
                                   i]
    tmp.min <- updated.single.mat[1, i]
    
    
    if (var(tmp$inhibition, na.rm = TRUE) == 0) {
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }
    #modified Jing 27/08. Original was fitting L.4 function by default, but it sohuld have been Ll.4. 
    if (nrow(tmp) == 1) {
      fitted.inhibition = response.mat[c(2:nrow(response.mat)), i]
    } else {
      tmp.model = tryCatch(    
        {tmp.model <- drm(inhibition ~ dose, data = tmp, fct = LL.4(fixed = c(NA, 
                                                                                                   tmp.min, Emax, NA)), na.action = na.omit)
    }, warning = function(w) {
      tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, 
                                                                          tmp.min, Emax, NA)), na.action = na.omit)
    }, error = function(e) {
      tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, 
                                                                          tmp.min, Emax, NA)), na.action = na.omit)
    }
      )
      fitted.inhibition <- suppressWarnings(fitted(tmp.model))
    }
#original
#        tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, 
#                                                                        tmp.min, Emax, NA)), na.action = na.omit)
   
    
   #original 
  #   tmp$fitted.inhibition <- suppressWarnings(fitted(tmp.model))
    
    #modified by jing 27/08
    tmp$fitted.inhibition = fitted.inhibition
     
    # if (tmp$fitted.inhibition[nrow(response.mat) - 1] < 0) 
    #  tmp$fitted.inhibition[nrow(response.mat) - 1] <- tmp.min
    updated.col.mat[c(2:nrow(response.mat)), i] <- tmp$fitted.inhibition
  }
  updated.row.mat <- updated.single.mat
  for (i in 2:nrow(response.mat)) {
    tmp <- as.data.frame(mat.or.vec(ncol(response.mat) - 
                                      1, 0))
    tmp$dose <- as.numeric(colnames(response.mat)[-1])
    tmp$inhibition <- response.mat[i, c(2:ncol(response.mat))]
    tmp.min <- updated.single.mat[i, 1]
    if (var(tmp$inhibition, na.rm = TRUE) == 0) {
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }
    
    
    if (nrow(tmp) == 1) {
      fitted.inhibition = response.mat[i, c(2:ncol(response.mat))]
    } else {
      tmp.model = tryCatch(    
        {tmp.model <- drm(inhibition ~ dose, data = tmp, fct = LL.4(fixed = c(NA, 
                                                                              tmp.min, Emax, NA)), na.action = na.omit)
        }, warning = function(w) {
          tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, 
                                                                              tmp.min, Emax, NA)), na.action = na.omit)
        }, error = function(e) {
          tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, 
                                                                              tmp.min, Emax, NA)), na.action = na.omit)
        }
      )
      fitted.inhibition <- suppressWarnings(fitted(tmp.model))
    }
    tmp$fitted.inhibition <- fitted.inhibition
    updated.row.mat[i, c(2:ncol(response.mat))] <- tmp$fitted.inhibition
    }
    #tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, 
    #                                                                    tmp.min, Emax, NA)), na.action = na.omit)
   
    
#     tmp$fitted.inhibition <- suppressWarnings(fitted(tmp.model))
 #   if (tmp$fitted.inhibition[ncol(response.mat) - 1] < 0) 
#      tmp$fitted.inhibition[ncol(response.mat) - 1] <- tmp.min
 #   updated.row.mat[i, c(2:ncol(response.mat))] <- tmp$fitted.inhibition
#  }
  fitted.mat <- (updated.col.mat + updated.row.mat)/2
  zip.mat <- updated.single.mat
  for (i in 2:nrow(updated.single.mat)) {
    for (j in 2:ncol((updated.single.mat))) {
      zip.mat[i, j] <- updated.single.mat[i, 1] + updated.single.mat[1, 
                                                                     j] - updated.single.mat[i, 1] * updated.single.mat[1, 
                                                                                                                        j]/100
    }
  }
  fitted.mat[1, 1] <- 0
  zip.mat[1, 1] <- 0
  fitted.mat <- apply(fitted.mat, c(1, 2), function(x) ifelse(x > 
                                                                100, 100, x))
  delta.mat <- (fitted.mat - zip.mat)
  delta.mat
}

ZIP <- function (response.mat,correction = T, Emin = 0, Emax = 100,nan.handle = c("LL4", "L4")) {
  # scores <- list()
  #  method <- "Loewe"
  out <- tryCatch(ZIP_modded(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
                             nan.handle = nan.handle), error = function(e) NA)
  #  {
  #  x$method <- method
  #  x$scores <- x$dose.response.mats
  #  x$scores[!is.na(x$scores)] <- NA
  #  x$scores <- replace(x$scores, !is.na(x$scores), NA)
  #  purrr::map(!is.na(x$scores), NA)
  #  x$scores <- lapply(x$scores, function(z) ifelse(!is.na(z), NA, z))
  #  return(x)
  #  })
  return(out)
}

BLISS_modded<- function (response.mat, correction = T, Emin = 0, Emax = 100, 
                         nan.handle = c("LL4", "L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD(response.mat, Emin = Emin, 
                                         Emax = Emax, nan.handle)$corrected.mat
  }
  drug1.response <- response.mat[, 1]
  drug2.response <- response.mat[1, ]
  ref.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      ref.mat[i, j] <- drug1.response[i] + drug2.response[j] - 
        drug1.response[i] * drug2.response[j]/100
    }
  }
  syn.mat <- response.mat - ref.mat
  syn.mat
}

Bliss <- function (response.mat,correction =T, Emin = 0, Emax = 100,nan.handle = c("LL4", "L4")) {
  # scores <- list()
  #  method <- "ZIP"
  out <- tryCatch(BLISS_modded(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
                               nan.handle = nan.handle), error = function(e) NA)
  #  {
  #  x$method <- method
  #  x$scores <- x$dose.response.mats
  #  x$scores[!is.na(x$scores)] <- NA
  #  x$scores <- replace(x$scores, !is.na(x$scores), NA)
  #  purrr::map(!is.na(x$scores), NA)
  #  x$scores <- lapply(x$scores, function(z) ifelse(!is.na(z), NA, z))
  #  return(x)
  #  })
  return(out)
}

HSA <- function (response.mat,correction = T, Emin = 0, Emax = 100,nan.handle = c("LL4", "L4")) {
  # scores <- list()
  #  method <- "ZIP"
  out <- tryCatch(HSA_modded(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
                             nan.handle = nan.handle), error = function(e) NA)
  #  {
  #  x$method <- method
  #  x$scores <- x$dose.response.mats
  #  x$scores[!is.na(x$scores)] <- NA
  #  x$scores <- replace(x$scores, !is.na(x$scores), NA)
  #  purrr::map(!is.na(x$scores), NA)
  #  x$scores <- lapply(x$scores, function(z) ifelse(!is.na(z), NA, z))
  #  return(x)
  #  })
  return(out)
}

HSA_modded <-function (response.mat, correction = T, Emin = 0, Emax = 100, 
                       nan.handle = c("LL4", "L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD(response.mat, Emin = Emin, 
                                         Emax = Emax, nan.handle)$corrected.mat
  }
  drug1.response <- response.mat[, 1]
  drug2.response <- response.mat[1, ]
  ref.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      ref.mat[i, j] <- ifelse(drug1.response[i] > drug2.response[j], 
                              drug1.response[i], drug2.response[j])
    }
  }
  syn.mat <- response.mat - ref.mat
  syn.mat
}

#original
#CalculateSynergy <- function (data, method = "ZIP", correction = TRUE, Emin = 0, 
#                              Emax = 100, nan.handle = c("LL4", "L4")) 
# modified by Jing 27/08
CalculateSynergy <- function (data, method = "ZIP", correction = TRUE, Emin = 0, Emax = NA, nan.handle = c("L4"))
{
  if (!is.list(data)) {
    stop("Input data is not a list format!")
  }
  if (!method %in% c("ZIP", "HSA", "Bliss", "Loewe")) {
    stop("The method parameter can only be one of the following: ZIP, HSA, Bliss and Loewe.")
  }
  dose.response.mats <- data$dose.response.mats
  num.pairs <- length(dose.response.mats)
  scores <- list()
  nan.handle <- match.arg(nan.handle)
  for (i in 1:num.pairs) {
    response.mat <- dose.response.mats[[i]]
    scores[[i]] <- switch(method, ZIP = ZIP(response.mat, 
                                            correction, Emin = Emin, Emax = Emax, nan.handle), 
                          HSA = HSA(response.mat, correction, Emin = Emin, 
                                    Emax = Emax, nan.handle), Bliss = Bliss(response.mat, 
                                                                            correction, Emin = Emin, Emax = Emax, nan.handle), 
                          Loewe = Loewe(response.mat, correction, Emin = Emin, 
                                        Emax = Emax, nan.handle))
  }
  data$scores <- scores
  data$method <- method
  return(data)
}

BaselineCorrectionSD <- function (response.mat, Emin = NA, Emax = NA, nan.handle = c("LL4", 
                                                             "L4")) 
{
  pm <- response.mat
  if (is.null(rownames(response.mat)) | is.null(colnames(response.mat))) {
    stop("Please provide drug contrations as row names and column names!")
  }
  nan.handle <- match.arg(nan.handle)
  single.fitted <- FittingSingleDrug(response.mat, c(NA, Emin, 
                                                     Emax, NA), nan.handle)
  baseline <- (min(as.numeric(single.fitted$drug.row.fitted)) + 
                 min(as.numeric(single.fitted$drug.col.fitted)))/2
  pm.corrected <- pm - ((100 - pm)/100 * baseline)
  output <- list(original.mat = pm, corrected.mat = pm.corrected)
  return(output)
}

FittingSingleDrug <- function (response.mat, fixed = c(NA, NA, NA, NA), nan.handle = c("LL4", 
                                                                                       "L4")) 
{
  r.num <- nrow(response.mat)
  c.num <- ncol(response.mat)
  #drug.col <- cbind(as.numeric(colnames(response.mat)[-1]), #original
  drug.col <- cbind(as.numeric(colnames(response.mat)), # modified by Jing 27/08
 #                   response.mat[1, 2:c.num]) #original
                      response.mat[1, 1:c.num]) # modified by Jing 27/08
 colnames(drug.col) <- c("conc", "effect")
  drug.col <- as.data.frame(apply(drug.col, 2, as.numeric))
 
  #modified JIng 27/08 
  if(nrow(drug.col) != 1) {
    if (var(drug.col$effect) == 0) {
    drug.col$effect[nrow(drug.col)] <- drug.col$effect[nrow(drug.col)] + 
      10^-10
    }}
 #original 
  #if (var(drug.col$effect) == 0) {
  #  drug.col$effect[nrow(drug.col)] <- drug.col$effect[nrow(drug.col)] + 
  #    10^-10
  #}
  
  
  
  nan.handle <- match.arg(nan.handle)
  drug.col.model <- tryCatch({
    drm(effect ~ conc, data = drug.col, fct = LL.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE))
  }, warning = function(w) {
    if (nan.handle == "L4") {
      drm(effect ~ conc, data = drug.col, fct = L.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE))
    }
    else {
      drm(effect ~ conc, data = drug.col, fct = LL.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE))
    }
  }, error = function(e) {
    drm(effect ~ conc, data = drug.col, fct = L.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE))
  })
  drug.col.fitted <- suppressWarnings(fitted(drug.col.model))
  
  #original
 # drug.row <- cbind(as.numeric(rownames(response.mat)[-1]), 
#                    response.mat[2:r.num, 1])

  #modified Jing 27/08
  drug.row <- cbind(as.numeric(rownames(response.mat)), 
                    response.mat[1:r.num, 1])
  
  colnames(drug.row) <- c("conc", "effect")
  drug.row <- as.data.frame(apply(drug.row, 2, as.numeric))
#original
 #   if (var(drug.row$effect) == 0) { 
#    drug.row$effect[nrow(drug.row)] <- drug.row$effect[nrow(drug.row)] + 
#      10^-10
#    }
  #modified by Jing 27/08
  if (nrow(drug.row) != 1) {
    if (var(drug.row$effect) == 0) { 
      drug.row$effect[nrow(drug.row)] <- drug.row$effect[nrow(drug.row)] + 
        10^-10
    }
  }
  
  drug.row.model <- tryCatch({
    drm(effect ~ conc, data = drug.row, fct = LL.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE))
  }, warning = function(w) {
    if (nan.handle == "L4") {
      drm(effect ~ conc, data = drug.row, fct = L.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE))
    }
    else {
      drm(effect ~ conc, data = drug.row, fct = LL.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE))
    }
  }, error = function(e) {
    drm(effect ~ conc, data = drug.row, fct = L.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE))
  })
  drug.row.fitted <- suppressWarnings(fitted(drug.row.model))
  return(list(drug.row.fitted = drug.row.fitted, drug.row.model = drug.row.model, 
              drug.col.model = drug.col.model, drug.col.fitted = drug.col.fitted))
}

ReshapeData <- function (data, data.type = "viability") 
{
  if (!all(c("BlockID", "DrugRow", "DrugCol", "Row", "Col", 
             "Response", "ConcRow", "ConcCol", "ConcRowUnit", "ConcColUnit") %in% 
           colnames(data))) 
    stop("The input data must contain the following columns: BlockID, DrugRow, DrugCol, Row, Col, Response,\n         ConcRow, ConcCol, ConcRowUnit, ConcColUnit")
  id.drug.comb <- unique(data$BlockID)
  dose.response.mats <- list()
  drug.pairs <- data.frame(drug.row = character(length(id.drug.comb)), 
                           drug.col = character(length(id.drug.comb)), concRUnit = character(length(id.drug.comb)), 
                           concCUnit = character(length(id.drug.comb)), blockIDs = numeric(length(id.drug.comb)), 
                           stringsAsFactors = FALSE)
  for (i in 1:length(id.drug.comb)) {
    tmp.mat <- data[which(data$BlockID == id.drug.comb[i]), 
                    ]
    if (data.type == "viability") {
      tmp.mat$Inhibition <- 100 - tmp.mat$Response
    }
    else {
      tmp.mat$Inhibition <- tmp.mat$Response
    }
    conc.col <- tmp.mat$ConcCol[which(tmp.mat$Row == 1)]
    conc.col <- conc.col[order(tmp.mat$Col[which(tmp.mat$Row == 
                                                   1)])]
    conc.row <- tmp.mat$ConcRow[which(tmp.mat$Col == 1)]
    conc.row <- conc.row[order(tmp.mat$Row[which(tmp.mat$Col == 
                                                   1)])]
    response.mat <- acast(tmp.mat, Row ~ Col, value.var = "Inhibition")
    colnames(response.mat) <- conc.col
    rownames(response.mat) <- conc.row
    if (which.max(conc.row) == 1 & which.max(conc.col) == 
        1) {
      response.mat <- t(apply(apply(response.mat, 2, rev), 
                              1, rev))
    }
    else if (which.max(conc.row) == length(conc.row) & which.max(conc.col) == 
             1) {
      response.mat <- t(apply(response.mat, 1, rev))
    }
    else if (which.max(conc.row) == 1 & which.max(conc.col) == 
             length(conc.col)) {
      response.mat <- apply(response.mat, 2, rev)
    }
    conc.runit <- unique(tmp.mat$ConcRowUnit)
    conc.cunit <- unique(tmp.mat$ConcColUnit)
    drug.row <- unique(tmp.mat$DrugRow)
    drug.col <- unique(tmp.mat$DrugCol)
    drug.pairs$drug.row[i] <- drug.row
    drug.pairs$drug.col[i] <- drug.col
    drug.pairs$concRUnit[i] <- conc.runit
    drug.pairs$concCUnit[i] <- conc.cunit
    dose.response.mats[[i]] <- response.mat
  }
  drug.pairs$blockIDs <- id.drug.comb
  return(list(dose.response.mats = dose.response.mats, drug.pairs = drug.pairs))
}

#####
#some error handling
# this throws the dim error
#  3053  750690 in cell line SK-MEL-28
test <- for_sorting %>%
  filter(CELLNAME == 'SK-MEL-28')
x <- test
y <- 256942
z <- 256439


y <- as.numeric(y)
z <- as.numeric(z)

x <- filter(x, ((NSC1 == y & NSC2 == z) | (NSC1 == y & is.na(NSC2)) | (NSC1 == z & is.na(NSC2)) ))

x_comb <- list()
x_comb <- filter(x, ( (NSC1 == y & NSC2 == z) ) )
x_comb <- x_comb %>%
  group_by(NSC1, NSC2, CONC1, CONC2) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()


aCONC11 <- unique(filter(x_comb, (!is.na(CONC2))  )$CONC1)
aNSC11 <- y
aCONC21 <- unique(filter(x_comb, (!is.na(CONC2))  )$CONC2)
aNSC21 <-z

x_one <- filter(x, (CONC1 %in% aCONC11)& NSC1 ==aNSC11 & is.na(CONC2))
x_one <- x_one %>%
  group_by(CONC1) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()
df <- data.frame()
df <- as.data.frame(cbind(order(aCONC11), aCONC11))
x_one <- inner_join(x_one, df, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC11"))
x_one <- x_one[!duplicated(x_one[c("CONC1")]),]

x_two <- filter(x, (CONC1 %in% aCONC21)& NSC1 ==aNSC21 & is.na(CONC2))
x_two <- x_two %>%
  group_by(CONC1) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()
df2 <- data.frame()
df2 <- as.data.frame(cbind(order(aCONC21), aCONC21))
x_two <- inner_join(x_two, df2, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC21"))
x_two <- x_two[!duplicated(x_two[c("CONC1")]),]

x_comb <- rbind(x_comb, x_one)
x_comb <- rbind(x_comb, x_two)

x <- x_comb

class_old <- list()
class_old <- sapply(x, class)
#  print(class_old)
NSC2_un <- vector(mode="numeric", length=1)
NSC2_un <- unique(x$NSC2[!is.na(x$NSC2)])
predicate <- function(x) {
  if (as.numeric(x[["CONCINDEX2"]]== 0) & as.numeric(x[["NSC1"]]) == NSC2_un) {
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
#  x <- x[,-c(1:16)] # to-do change that
x <- x[-match(c("COMBODRUGSEQ","SCREENER","STUDY", "TESTDATE", "PLATE", "PANELNBR", "CELLNBR", "PREFIX1", "NSC1", "SAMPLE1", "CONCINDEX1", "CONC1", "CONCUNIT1", "PREFIX2", "NSC2", "SAMPLE2", "CONCINDEX2", "CONC2", "CONCUNIT2", "PERCENTGROWTH", "PERCENTGROWTHNOTZ","TESTVALUE", "CONTROLVALUE", "TZVALUE", "EXPECTEDGROWTH", "SCORE", "VALID", "PANEL", "CELLNAME", "mean_single"), names(x))]
return(x)


#####

do.call(rbind, datalist_full) -> trying
do.call(rbind, trying) -> trying

#######################################
function (data, data.type = "viability") 
{
  if (!all(c("BlockID", "DrugRow", "DrugCol", "Row", "Col", 
             "Response", "ConcRow", "ConcCol", "ConcRowUnit", "ConcColUnit") %in% 
           colnames(data))) 
    stop("The input data must contain the following columns: BlockID, DrugRow, DrugCol, Row, Col, Response,\n         ConcRow, ConcCol, ConcRowUnit, ConcColUnit")
  id.drug.comb <- unique(data$BlockID)
  dose.response.mats <- list()
  drug.pairs <- data.frame(drug.row = character(length(id.drug.comb)), 
                           drug.col = character(length(id.drug.comb)), concRUnit = character(length(id.drug.comb)), 
                           concCUnit = character(length(id.drug.comb)), blockIDs = numeric(length(id.drug.comb)), 
                           stringsAsFactors = FALSE)
  for (i in 1:length(id.drug.comb)) {
    tmp.mat <- data[which(data$BlockID == id.drug.comb[i]), 
                    ]
    if (data.type == "viability") {
      tmp.mat$Inhibition <- 100 - tmp.mat$Response
    }
    else {
      tmp.mat$Inhibition <- tmp.mat$Response
    }
    conc.col <- tmp.mat$ConcCol[which(tmp.mat$Row == 1)]
    conc.col <- conc.col[order(tmp.mat$Col[which(tmp.mat$Row == 
                                                   1)])]
    conc.row <- tmp.mat$ConcRow[which(tmp.mat$Col == 1)]
    conc.row <- conc.row[order(tmp.mat$Row[which(tmp.mat$Col == 
                                                   1)])]
    response.mat <- acast(tmp.mat, Row ~ Col, value.var = "Inhibition")
    colnames(response.mat) <- conc.col
    rownames(response.mat) <- conc.row
    if (which.max(conc.row) == 1 & which.max(conc.col) == 
        1) {
      response.mat <- t(apply(apply(response.mat, 2, rev), 
                              1, rev))
    }
    else if (which.max(conc.row) == length(conc.row) & which.max(conc.col) == 
             1) {
      response.mat <- t(apply(response.mat, 1, rev))
    }
    else if (which.max(conc.row) == 1 & which.max(conc.col) == 
             length(conc.col)) {
      response.mat <- apply(response.mat, 2, rev)
    }
    conc.runit <- unique(tmp.mat$ConcRowUnit)
    conc.cunit <- unique(tmp.mat$ConcColUnit)
    drug.row <- unique(tmp.mat$DrugRow)
    drug.col <- unique(tmp.mat$DrugCol)
    drug.pairs$drug.row[i] <- drug.row
    drug.pairs$drug.col[i] <- drug.col
    drug.pairs$concRUnit[i] <- conc.runit
    drug.pairs$concCUnit[i] <- conc.cunit
    dose.response.mats[[i]] <- response.mat
  }
  drug.pairs$blockIDs <- id.drug.comb
  return(list(dose.response.mats = dose.response.mats, drug.pairs = drug.pairs))
}



#######################################  
testing <- for_sorting %>%
  filter(CELLNAME == "HCT-116")

y = 681239
z =  606869
y = 3053
z = 750690
testing <- filter(big, ((DrugCol == y & DrugRow == z) | (DrugCol == y & is.na(DrugRow)) | (DrugRow == z & is.na(DrugCol))  ))

testing1 <- filter(testing, ( (NSC1 == y & NSC2 == z) ) )
testing1 <- testing1 %>%
  group_by(NSC1, NSC2, CONC1, CONC2) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()

aCONC11 <- unique(filter(testing1, (!is.na(CONC2))  )$CONC1)
aNSC11 <- y #unique(filter(x, (!is.na(CONC2))  )$NSC1)
aCONC21 <- unique(filter(testing1, (!is.na(CONC2))  )$CONC2)
aNSC21 <-z

testing2 <- filter(testing, (CONC1 %in% aCONC11)& NSC1 ==aNSC11 & is.na(CONC2))
testing2 <- testing2 %>%
  group_by(CONC1) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()
df <- as.data.frame(cbind(order(aCONC11), aCONC11))
inner_join(testing2, df, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC11")) -> testing6
testing6 <- testing6[!duplicated(testing6[c("CONC1")]),]

testing3 <- filter(testing, (CONC1 %in% aCONC21)& NSC1 ==aNSC21 & is.na(CONC2))
testing3 <- testing3 %>%
  group_by(CONC1) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()
df2 <- as.data.frame(cbind(order(aCONC21), aCONC21))
inner_join(testing3, df2, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC21")) -> testing7
testing7 <- testing7[!duplicated(testing7[c("CONC1")]),]

testing1 <- rbind(testing1, testing6)
testing1 <- rbind(testing1, testing7)


rowInsert(testing1, 'test') -> testing7


testing <- testing %>%
  group_by(CONCINDEX1,CONC1, CONC2, NSC1, NSC2) %>%
  mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
  arrange()

aCONC1 <- unique(filter(testing, (!is.na(CONC2))  )$CONC1)
aNSC1 <- y #unique(filter(x, (!is.na(CONC2))  )$NSC1)
aCONC2 <- unique(filter(testing, (!is.na(CONC2))  )$CONC2)
aNSC2 <-z

testing <- filter(testing, ( (!is.na(NSC1) & !is.na(NSC2)) | ((CONC1 %in% aCONC1)& NSC1 ==aNSC1) | ((CONC1 %in% aCONC2)&NSC1==aNSC2) )  )

testing <- testing[!duplicated(testing[c("NSC1","NSC2" ,"CONC1", "CONC2", "mean_single")]),]

screenID <- prettyR::Mode(testing$SCREENER)

testing <- filter(testing,SCREENER ==  screenID)

clorafarbine 759857
bortezomib 681239