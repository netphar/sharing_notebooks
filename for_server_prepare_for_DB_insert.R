library("tidyverse")
library('drc')
library('reshape2')
library('nleqslv')

#for debug on laptop
#setwd('/Users/zagidull/Documents/fimm_files/synergy_calc_august/almanac')
#input <- readRDS(file = '2208_first_100_combos')

#for debug and production (lol) on server
setwd('/home/bulat/NCI/almanac')
reshaped.Bliss <- readRDS(file = 'reshaped.Bliss')
reshaped.ZIP <- readRDS(file = 'reshaped.ZIP')
reshaped.HSA <- readRDS(file = 'reshaped.HSA')
reshaped.Loewe <- readRDS(file = 'reshaped.Loewe')



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
first.two$dose.response.mats <- reshaped.ZIP$dose.response.mats
first.two$drug.pairs <- reshaped.ZIP$drug.pairs
first.two$scores <- reshaped.ZIP$scores
first.two$method <- reshaped.ZIP$method

#populate for Bliss
first.two.Bliss$dose.response.mats <- reshaped.Bliss$dose.response.mats
first.two.Bliss$drug.pairs <- reshaped.Bliss$drug.pairs
first.two.Bliss$scores <- reshaped.Bliss$scores
first.two.Bliss$method <- reshaped.Bliss$method

#populate for HSA
first.two.HSA$dose.response.mats <- reshaped.HSA$dose.response.mats
first.two.HSA$drug.pairs <- reshaped.HSA$drug.pairs
first.two.HSA$scores <- reshaped.HSA$scores
first.two.HSA$method <- reshaped.HSA$method

#populate for Loewe
first.two.Loewe$dose.response.mats <- reshaped.Loewe$dose.response.mats
first.two.Loewe$drug.pairs <- reshaped.Loewe$drug.pairs
first.two.Loewe$scores <- reshaped.Loewe$scores
first.two.Loewe$method <- reshaped.Loewe$method

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
saveRDS(object = big_data, file = 'big_data_with_Loewe')
saveRDS(object = datalist, file = 'datalist_with_Loewe')
