# this script is used to run CalculateSynergy on fimm server
# version of CalculateSynergy is current as of 28.08
# input is R list of lists, result synergyfinder::ReshapeData()
# NB: there is error handling added in form of tryCatch for CalculateSynergy. As well as in downstream functions
# NB1: global options(warn = -1) and options(warn = 0) are set throughout the script.
# input with full list of valid dose.response.mats and corresponding drugs pairs is 2308_reshaped_finished_4950688

#housekeeping
rm(list=ls(all=TRUE))

#libs
library("tidyverse")
library('drc')
library('reshape2')
library('nleqslv')

#for debug on laptop
#setwd('/Users/zagidull/Documents/fimm_files/synergy_calc_august/almanac')
#input <- readRDS(file = '2208_first_100_combos')

#for debug and production (lol) on server
#setwd('/home/bulat/NCI/almanac')
input <- readRDS(file = '/Users/zagidull/Documents/fimm_files/synergy_calc_august/almanac/2308_reshaped_finished_4950688')

#input$dose.response.mats <- input$dose.response.mats[c(1:3)]
#input$drug.pairs <- input$drug.pairs[c(1:3),]


#change concetrations from M to uM
input$dose.response.mats <- lapply(input$dose.response.mats, function(x) {
  colnames(x) <- as.numeric(colnames(x))*(10^6)
  rownames(x) <- as.numeric(rownames(x))*(10^6)
  return(x)
})
input$drug.pairs$concRUnit <- input$drug.pairs$concCUnit <- 'uM'
######end of mod from 2/09

HSA <- function (response.mat,correction = T, Emin = 0, Emax = NA,nan.handle = c("L4")) {
  # scores <- list()
  #  method <- "ZIP"
  out <- tryCatch(HSA2(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
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

HSA2 = function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
                 nan.handle = c("L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD2(response.mat, Emin = Emin, 
                                          Emax = Emax, nan.handle = nan.handle)$corrected.mat
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

#jing's baseline correction
BaselineCorrectionSD2 = function (response.mat, Emin = NA, Emax = NA, nan.handle = c("L4")) 
{
  pm <- response.mat
  if (is.null(rownames(response.mat)) | is.null(colnames(response.mat))) {
    stop("Please provide drug contrations as row names and column names!")
  }
  nan.handle <- match.arg(nan.handle)
  single.fitted <- FittingSingleDrug2(response.mat, c(NA, Emin, 
                                                      Emax, NA), nan.handle)
  baseline <- (min(as.numeric(single.fitted$drug.row.fitted)) + 
                 min(as.numeric(single.fitted$drug.col.fitted)))/2
  pm.corrected <- pm - ((100 - pm)/100 * baseline)
  output <- list(original.mat = pm, corrected.mat = pm.corrected)
  return(output)
}

#jing's fitting single
FittingSingleDrug2 = function (response.mat, fixed = c(NA, NA, NA, NA), nan.handle = c("L4")) 
{
  r.num <- nrow(response.mat)
  c.num <- ncol(response.mat)
  drug.col <- cbind(as.numeric(colnames(response.mat)), 
                    response.mat[1, 1:c.num])
  drug.col <- data.frame(drug.col)
  colnames(drug.col) <- c("conc", "effect")
  drug.col$conc = as.numeric(drug.col$conc)
  drug.col$effect = as.numeric(drug.col$effect)
  
  if(nrow(drug.col)!=1){
    if (var(drug.col$effect) == 0) {
      drug.col$effect[nrow(drug.col)] <- drug.col$effect[nrow(drug.col)] + 
        10^-10
    }
  }
  
  nan.handle <- match.arg(nan.handle)
  drug.col.model <- tryCatch({
    drm(effect ~ conc, data = drug.col, fct = LL.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
  }, warning = function(w) {
    if (nan.handle == "L4") {
      drm(effect ~ conc, data = drug.col, fct = L.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
    }
    else {
      drm(effect ~ conc, data = drug.col, fct = LL.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
    }
  }, error = function(e) {
    drm(effect ~ conc, data = drug.col, fct = L.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
  })
  drug.col.fitted <- suppressWarnings(fitted(drug.col.model))
  
  
  drug.row <- cbind(as.numeric(rownames(response.mat)), 
                    response.mat[1:r.num, 1])
  drug.row <- data.frame(drug.row)
  colnames(drug.row) <- c("conc", "effect")
  
  drug.row$conc = as.numeric(drug.row$conc)
  drug.row$effect = as.numeric(drug.row$effect)
  
  if (nrow(drug.row)!=1){
    if (var(drug.row$effect) == 0) {
      drug.row$effect[nrow(drug.row)] <- drug.row$effect[nrow(drug.row)] + 
        10^-10
    }
  }
  drug.row.model <- tryCatch({
    drm(effect ~ conc, data = drug.row, fct = LL.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
  }, warning = function(w) {
    if (nan.handle == "L4") {
      drm(effect ~ conc, data = drug.row, fct = L.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
    }
    else {
      drm(effect ~ conc, data = drug.row, fct = LL.4(fixed = fixed), 
          na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
    }
  }, error = function(e) {
    drm(effect ~ conc, data = drug.row, fct = L.4(fixed = fixed), 
        na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
  })
  drug.row.fitted <- suppressWarnings(fitted(drug.row.model))
  
  res = list(drug.row.fitted = drug.row.fitted, drug.row.model = drug.row.model, 
             drug.col.model = drug.col.model, drug.col.fitted = drug.col.fitted)
  return(res)
}

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
  pb <- txtProgressBar(min = 0, max = num.pairs, style = 3)
  
  for (i in 1:num.pairs) {
    setTxtProgressBar(pb, i)
    
    response.mat <- dose.response.mats[[i]]
    scores[[i]] <- switch(method, ZIP = ZIP(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle), 
                          HSA = HSA(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle), 
                          Bliss = Bliss(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle), 
                          Loewe = Loewe(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle))
  }
  data$scores <- scores
  data$method <- method
  close(pb)
  
  return(data)
}

CalculateSynergy(input, method = 'HSA', correction = T, Emin = 0, Emax = NA) -> input.HSA
saveRDS(object = input.HSA, file = '0309_reshaped.HSA')


#CalculateSynergy(input, method = 'Bliss', correction = T, Emin = 0) -> input.Bliss
#saveRDS(object = input.Bliss, file = '3008_reshaped.Bliss')

