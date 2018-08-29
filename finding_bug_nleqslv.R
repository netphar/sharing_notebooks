library('drc')
library('reshape2')
library('nleqslv')

#Jing testing 29/08
res = readRDS("2308_reshaped_finished_4950688")
res$dose.response.mats[[1]]
res$drug.pairs[1,]
try <- list()
try$dose.response.mats <- res$dose.response.mats[1]
try$drug.pairs <-  res$drug.pairs[1,]
try.Bliss <- CalculateSynergy(try, method = 'Bliss', correction = T)
try.Loewe <- CalculateSynergy(try, method = 'Loewe', correction = T)
try.ZIP <- CalculateSynergy(try, method = 'ZIP', correction = T)
try.HSA <- CalculateSynergy(try, method = 'HSA', correction = T)

#jing's single fitted


#jing's loewe
Loewe2 = function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
                   nan.handle = c("L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD2(response.mat, Emin = Emin, 
                                          Emax = Emax, nan.handle)$corrected.mat
  }
  single.fit <- FittingSingleDrug2(response.mat,fixed = c(NA, NA, NA, NA), nan.handle)
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
BaselineCorrectionSD2 = function (response.mat, Emin = NA, Emax = NA, nan.handle = c("LL4", 
                                                                                     "L4")) 
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


#reshaping finished results
datalist = list()
mylist <- list()
mylist.Bliss <- list()
mylist.HSA <- list()
mylist.Loewe <- list()

first.two <- list()
first.two.Bliss <- list()
first.two.HSA <- list()
first.two.Loewe <- list()

#populate for ZIP
first.two$dose.response.mats <- try.ZIP$dose.response.mats
first.two$drug.pairs <- try.ZIP$drug.pairs
first.two$scores <- try.ZIP$scores
first.two$method <- try.ZIP$method

#populate for Bliss
first.two.Bliss$dose.response.mats <- try.Bliss$dose.response.mats
first.two.Bliss$drug.pairs <- try.Bliss$drug.pairs
first.two.Bliss$scores <- try.Bliss$scores
first.two.Bliss$method <- try.Bliss$method

#populate for HSA
first.two.HSA$dose.response.mats <- try.HSA$dose.response.mats
first.two.HSA$drug.pairs <- try.HSA$drug.pairs
first.two.HSA$scores <- try.HSA$scores
first.two.HSA$method <- try.HSA$method

#populate for Loewe
first.two.Loewe$dose.response.mats <- try.Loewe$dose.response.mats
first.two.Loewe$drug.pairs <- try.Loewe$drug.pairs
first.two.Loewe$scores <- try.Loewe$scores
first.two.Loewe$method <- try.Loewe$method

for (i in 1:nrow(first.two$drug.pairs)) {
  mylist[i] <- first.two$drug.pairs[i,]$blockIDs
}

for (i in 1:nrow(first.two.Bliss$drug.pairs)) {
  mylist.Bliss[i] <- first.two.Bliss$drug.pairs[i,]$blockIDs
}
for (i in 1:nrow(first.two.HSA$drug.pairs)) {
  mylist.HSA[i] <- first.two.HSA$drug.pairs[i,]$blockIDs
}

for (i in 1:nrow(first.two.Loewe$drug.pairs)) {
  mylist.Loewe[i] <- first.two.Loewe$drug.pairs[i,]$blockIDs
}

pb <- txtProgressBar(min = 0, max = nrow(first.two$drug.pairs), style = 3)

for (i in 1:nrow(first.two$drug.pairs))
{
  setTxtProgressBar(pb, i)
  a <- list()
  b <- list()
  c <- list()
  d <- list()
  
  
  a$dose.response.mats <- first.two$dose.response.mats[i]
  a$scores <- first.two$scores[i]
  a$method <- first.two$method
  a$drug.pairs <- first.two$drug.pairs[i,]
  
  if (is.na(a$scores)) {
    cols <- colnames(a$dose.response.mats[[which(is.na(a$scores))]])
    rows <- rownames(a$dose.response.mats[[which(is.na(a$scores))]])
    numrows <- dim(a$dose.response.mats[[which(is.na(a$scores))]])[1]
    numcols <- dim(a$dose.response.mats[[which(is.na(a$scores))]])[2]
    a$scores[[which(is.na(a$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  b$dose.response.mats <- first.two.Bliss$dose.response.mats[i]
  b$scores <- first.two.Bliss$scores[i]
  b$method <- first.two.Bliss$method
  b$drug.pairs <- first.two.Bliss$drug.pairs[i,]
  
  if (is.na(b$scores)) {
    cols <- colnames(b$dose.response.mats[[which(is.na(b$scores))]])
    rows <- rownames(b$dose.response.mats[[which(is.na(b$scores))]])
    numrows <- dim(b$dose.response.mats[[which(is.na(b$scores))]])[1]
    numcols <- dim(c$dose.response.mats[[which(is.na(b$scores))]])[2]
    b$scores[[which(is.na(b$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  c$dose.response.mats <- first.two.Loewe$dose.response.mats[i]
  c$scores <- first.two.Loewe$scores[i]
  c$method <- first.two.Loewe$method
  c$drug.pairs <- first.two.Loewe$drug.pairs[i,]
  
  if (is.na(c$scores)) {
    
    cols <- colnames(c$dose.response.mats[[which(is.na(c$scores))]])
    rows <- rownames(c$dose.response.mats[[which(is.na(c$scores))]])
    numrows <- dim(c$dose.response.mats[[which(is.na(c$scores))]])[1]
    numcols <- dim(c$dose.response.mats[[which(is.na(c$scores))]])[2]
    c$scores[[which(is.na(c$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  d$dose.response.mats <- first.two.HSA$dose.response.mats[i]
  d$scores <- first.two.HSA$scores[i]
  d$method <- first.two.HSA$method
  d$drug.pairs <- first.two.HSA$drug.pairs[i,]
  
  if (is.na(d$scores)) {
    cols <- colnames(d$dose.response.mats[[which(is.na(d$scores))]])
    rows <- rownames(d$dose.response.mats[[which(is.na(d$scores))]])
    numrows <- dim(d$dose.response.mats[[which(is.na(d$scores))]])[1]
    numcols <- dim(d$dose.response.mats[[which(is.na(d$scores))]])[2]
    d$scores[[which(is.na(d$scores))]] <- matrix(,numrows,numcols,dimnames = list(rows, cols))
    
  }
  
  x.ZIP <- melt(a$dose.response.mats)
  y.ZIP <- melt(a$scores)
  xy.ZIP <- merge(x.ZIP,y.ZIP, by = c('Var1', 'Var2', 'L1'))
  
  x.Bliss <- melt(b$dose.response.mats)
  y.Bliss <- melt(b$scores)
  xy.Bliss <- merge(x.Bliss,y.Bliss, by = c('Var1', 'Var2', 'L1'))
  
  x.Loewe <- melt(c$dose.response.mats)
  y.Loewe <- melt(c$scores)
  xy.Loewe <- merge(x.Loewe,y.Loewe, by = c('Var1', 'Var2', 'L1'))
  
  x.HSA <- melt(d$dose.response.mats)
  y.HSA <- melt(d$scores)
  xy.HSA <- merge(x.HSA,y.HSA, by = c('Var1', 'Var2', 'L1'))
  
  colnames(xy.ZIP)[which(names(xy.ZIP) == "Var1")] <- "ConcR"
  colnames(xy.ZIP)[which(names(xy.ZIP) == "Var2")] <- "ConcC"
  colnames(xy.ZIP)[which(names(xy.ZIP) == "L1")] <- "blockIDs"
  xy.ZIP$blockIDs <- mylist[[i]]
  colnames(xy.ZIP)[which(names(xy.ZIP) == "value.x")] <- "Response_inhibition"
  synergy.type.ZIP <- paste('Synergy', first.two$method, sep = '_')
  colnames(xy.ZIP)[which(names(xy.ZIP) == "value.y")] <- synergy.type.ZIP
  
  colnames(xy.Bliss)[which(names(xy.Bliss) == "Var1")] <- "ConcR"
  colnames(xy.Bliss)[which(names(xy.Bliss) == "Var2")] <- "ConcC"
  colnames(xy.Bliss)[which(names(xy.Bliss) == "L1")] <- "blockIDs"
  xy.Bliss$blockIDs <- mylist.Bliss[[i]]
  colnames(xy.Bliss)[which(names(xy.Bliss) == "value.x")] <- "Response_inhibition"
  synergy.type.Bliss <- paste('Synergy', first.two.Bliss$method, sep = '_')
  colnames(xy.Bliss)[which(names(xy.Bliss) == "value.y")] <- synergy.type.Bliss
  
  colnames(xy.Loewe)[which(names(xy.Loewe) == "Var1")] <- "ConcR"
  colnames(xy.Loewe)[which(names(xy.Loewe) == "Var2")] <- "ConcC"
  colnames(xy.Loewe)[which(names(xy.Loewe) == "L1")] <- "blockIDs"
  xy.Loewe$blockIDs <- mylist.Loewe[[i]]
  colnames(xy.Loewe)[which(names(xy.Loewe) == "value.x")] <- "Response_inhibition"
  synergy.type.Loewe <- paste('Synergy', first.two.Loewe$method, sep = '_')
  colnames(xy.Loewe)[which(names(xy.Loewe) == "value.y")] <- synergy.type.Loewe
  
  colnames(xy.HSA)[which(names(xy.HSA) == "Var1")] <- "ConcR"
  colnames(xy.HSA)[which(names(xy.HSA) == "Var2")] <- "ConcC"
  colnames(xy.HSA)[which(names(xy.HSA) == "L1")] <- "blockIDs"
  xy.HSA$blockIDs <- mylist.HSA[[i]]
  colnames(xy.HSA)[which(names(xy.HSA) == "value.x")] <- "Response_inhibition"
  synergy.type.HSA <- paste('Synergy', first.two.HSA$method, sep = '_')
  colnames(xy.HSA)[which(names(xy.HSA) == "value.y")] <- synergy.type.HSA
  
  
  first.two.no.celllines <- merge(xy.ZIP, a$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  first.two.no.celllines.Bliss <- merge(xy.Bliss, b$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  first.two.no.celllines.Loewe <- merge(xy.Loewe, c$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  first.two.no.celllines.HSA <- merge(xy.HSA, d$drug.pairs, by.x=c('blockIDs'), by.y = c('blockIDs'), all.x = T)
  
  first.two.no.celllines.ZIP.Bliss <- merge(first.two.no.celllines, first.two.no.celllines.Bliss, 
                                            by = c('blockIDs','ConcR','ConcC','Response_inhibition',
                                                   'drug.row', 'drug.col','concRUnit','concCUnit'))
  first.two.no.celllines.ZIP.Bliss.HSA <- merge(first.two.no.celllines.ZIP.Bliss, first.two.no.celllines.HSA, 
                                                by = c('blockIDs','ConcR','ConcC','Response_inhibition',
                                                       'drug.row', 'drug.col','concRUnit','concCUnit'))
  first.two.no.celllines.ZIP.Bliss.HSA.Loewe <- merge(first.two.no.celllines.ZIP.Bliss.HSA, first.two.no.celllines.Loewe, 
                                                      by = c('blockIDs','ConcR','ConcC','Response_inhibition',
                                                             'drug.row', 'drug.col','concRUnit','concCUnit'))
  
  #  first.two.no.celllines.ZIP.Bliss.Loewe <- merge(first.two.no.celllines.ZIP.Bliss, first.two.no.celllines.Loewe, 
  #                                          by = c('blockIDs','ConcR','ConcC','Response_inhibition',
  #                                                 'drug.row', 'drug.col','concRUnit','concCUnit'))
  #  first.two.no.celllines.ZIP.Bliss.Loewe.HSA <- merge(first.two.no.celllines.ZIP.Bliss.Loewe, first.two.no.celllines.HSA, 
  #                                          by = c('blockIDs','ConcR','ConcC','Response_inhibition',
  #                                                 'drug.row', 'drug.col','concRUnit','concCUnit'))
  
  
  
  #  datalist[[i]] <- first.two.no.celllines.ZIP.Bliss.Loewe.HSA
  datalist[[i]] <- first.two.no.celllines.ZIP.Bliss.HSA.Loewe
  
}
close(pb)
big_data = do.call(rbind, datalist)





#to get data
tmp3 <- list()
readRDS('2308_reshaped_finished_4950688') -> reshaped.finished
from_reshaped <- reshaped.finished$drug.pairs[which(reshaped.finished[["drug.pairs"]]$blockIDs=='156:NCI-H322M'),]
reshaped.finished$dose.response.mats[4192] -> tmp3$dose.response.mats

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
Loewe_modded <- function (response.mat, correction = TRUE, Emin = 0, Emax = NA, 
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

ZIP_modded <- function (response.mat, correction = T, Emin = 0, Emax = NA, 
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

ZIP <- function (response.mat,correction = T, Emin = 0, Emax = NA,nan.handle = c("LL4", "L4")) {
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

Bliss <- function (response.mat,correction =T, Emin = 0, Emax = NA,nan.handle = c("LL4", "L4")) {
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

HSA <- function (response.mat,correction = T, Emin = 0, Emax = NA,nan.handle = c("LL4", "L4")) {
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

HSA_modded <-function (response.mat, correction = T, Emin = NA, Emax = NA, 
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