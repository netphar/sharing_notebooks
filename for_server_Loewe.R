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
setwd('/Users/zagidull/Documents/fimm_files/synergy_calc_august/almanac')
#input <- readRDS(file = '2208_first_100_combos')

#for debug and production (lol) on server
#setwd('/home/bulat/NCI/almanac')
input <- readRDS(file = '2308_reshaped_finished_4950688')
#input$dose.response.mats <- input$dose.response.mats[1]
#input$drug.pairs <- input$drug.pairs[1,]

input$dose.response.mats <- input$dose.response.mats[c(1:3)]
input$drug.pairs <- input$drug.pairs[c(1:3),]


#change concetrations from M to uM
input$dose.response.mats <- lapply(input$dose.response.mats, function(x) {
  colnames(x) <- as.numeric(colnames(x))*(10^6)
  rownames(x) <- as.numeric(rownames(x))*(10^6)
  return(x)
})
input$drug.pairs$concRUnit <- input$drug.pairs$concCUnit <- 'uM'
######end of mod from 2/09

Loewe <- function (response.mat,correction = T, Emin = 0, Emax = NA,nan.handle = c("L4")) {
  # scores <- list()
  #  method <- "Loewe"
  out <- tryCatch(Loewe2(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
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

Loewe2 = function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
                   nan.handle = c("L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD2(response.mat, Emin = Emin, 
                                          Emax = Emax, nan.handle = nan.handle)$corrected.mat
  }
  single.fit <- FittingSingleDrug2(response.mat,fixed = c(NA, NA, NA, NA), nan.handle = nan.handle)
  #  single.fit <- 
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

CalculateSynergy(input, method = 'Loewe', correction = T, Emin = 0, Emax = NA) -> input.Loewe
saveRDS(object = input.Loewe, file = '0309_reshaped.Loewe')


#CalculateSynergy(input, method = 'ZIP', correction = T, Emin = 0) -> input.ZIP
#saveRDS(object = input.ZIP, file = '3008_reshaped.ZIP')

#CalculateSynergy(input, method = 'HSA', correction = T, Emin = 0) -> input.HSA
#saveRDS(object = input.HSA, file = '3008_reshaped.HSA')


#CalculateSynergy(input, method = 'Bliss', correction = T, Emin = 0) -> input.Bliss
#saveRDS(object = input.Bliss, file = '3008_reshaped.Bliss')

