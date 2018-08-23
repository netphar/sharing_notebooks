library(dplyr)
library(synergyfinder)

#NB figure out how to sort out nesteded lists in place
# or use makeshift dictionary
# or just simple max() and extract the top n hits

#creates a list of key:value pairs. Where is key is WellID and value is Condition
mylist <- list()
for (i in seq_along(keys)) {
  mylist[keys[i]] <- values[i]
}

#this adds median scores sublist to the calculate synergy results 

CalculateSynergy(testing, method = "Bliss") -> testing.Bliss
CalculateSynergy(testing, method = "Loewe") -> testing.Loewe
CalculateSynergy(testing, method = "HSA") -> testing.HSA
CalculateSynergy(testing) -> testing.ZIP

CalculateMedianSynergy <- function(x) {
  x$median <- list()
  for (i in 1:nrow(x$drug.pairs)) {
  ifelse(is.na(x$scores[[i]]),
         x$median[[i]] <- NA, 
         x$median[[i]] <- round(mean(x$scores[[i]]), digits = 3) )
  }
  return(x)
}

CalculateMedianSynergy(testing.Bliss) -> testing.Bliss.with.medians
CalculateMedianSynergy(testing.Loewe) -> testing.Loewe.with.medians
CalculateMedianSynergy(testing.HSA) -> testing.HSA.with.medians
CalculateMedianSynergy(testing.ZIP) -> testing.ZIP.with.medians

#get top sorted median synergy scores
##sorts descending and creates a matrix with drug a, drug b, blockID, mean syn score
sorted.syn.scores.Bliss <- testing.Bliss.with.medians$median[order(testing.Bliss.with.medians$median, decreasing = T),]
top3 <- sorted.syn.scores[(1:3),]

#####
#end of working piece
####

#testing


    onego$median[[i]] <- round(mean(onego$scores[[i]]), digits = 3)
}
}
    
  current.score <- round(mean(onego$scores[[i]]), digits = 3)
  
  
  
  current.drugcombo <- workingcopy$drug.pairs[i,][c(1,2,5)]
  current.output <- cbind(current.drugcombo,current.score)
  mean.synergy.scores[[i]] <- current.output
  if (i == nrow(workingcopy$drug.pairs)) {
    output <- do.call(rbind, mean.synergy.scores)
    return(output)
  }
}


#creating a copy for current work
testing.ZIP -> workingcopy

#cleaning up NA scores
which(is.na(testing.ZIP$scores)) #get the id's of na scores

#positional list of na synergy scores
which(is.na(testing.ZIP$scores)) -> list.nas

#removing na scores 
workingcopy$scores[-c(which(is.na(testing.ZIP$scores)))] -> workingcopy$scores

#removing na matrices
workingcopy$dose.response.mats[-c(list.nas)] -> workingcopy$dose.response.mats

#removing na drug paires
workingcopy$drug.pairs[-c(list.nas),] -> workingcopy$drug.pairs

#trying drawing 3d response surfaces. NB: do not save files
synergyfinder::PlotSynergy(workingcopy,type = '3D', save.file = F)

#find out top 3 drug pairs with highest max ZIP synergy
mean.synergy.scores <- list()
for (i in 1:nrow(workingcopy$drug.pairs)) {
  current.score <- round(mean(workingcopy$scores[[i]]), digits = 3)
  current.drugcombo <- workingcopy$drug.pairs[i,][c(1,2,5)]
  current.output <- cbind(current.drugcombo,current.score)
  mean.synergy.scores[[i]] <- current.output
  if (i == nrow(workingcopy$drug.pairs)) {
   output <- do.call(rbind, mean.synergy.scores)
   return(output)
  }
}
##sorts descending and creates a matrix with drug a, drug b, blockID, mean syn score
sorted.syn.scores <- output[order(output$current.score, decreasing = T),]
top3 <- sorted.syn.scores[(1:3),]
##extract the blockID from there
top3.blockIDs <- top3[,3]

#draw plots for top 3 scores
for (blockIDs in top3.blockIDs) {
  
}


#figuring out top scores. Synergyfinder approach. NB: they just use mean(scores) rounded down to 3 sig figures
##setup for running the seeb elow bit
data <- workingcopy
scores <- data$scores
drug.pairs <- data$drug.pairs
num.pairs <- 1:length(scores)
#plots <- list()
#if (!is.null(pair.index)) {
  num.pairs <- pair.index
}
for (i in num.pairs) {
  scores.dose <- scores[[i]]
  drug.col <- drug.pairs$drug.col[i]
  drug.row <- drug.pairs$drug.row[i]
  scores.tmp <- scores.dose
  if (!is.null(col.range)) {
    if (col.range[1] == 1) {
      scores.tmp <- scores.tmp[, (col.range[1] + 1):col.range[2]]
    }
    else {
      scores.tmp <- scores.tmp[, col.range[1]:col.range[2]]
    }
  }
  else {
    scores.tmp <- scores.tmp[, -1]
  }
  if (!is.null(row.range)) {
    if (row.range[1] == 1) {
      scores.tmp <- scores.tmp[(row.range[1] + 1):row.range[2], 
                               ]
    }
    else {
      scores.tmp <- scores.tmp[row.range[1]:row.range[2], 
                               ]
    }
  }
  else {
    scores.tmp <- scores.tmp[-1, ]
  }
  summary.score <- round(mean(scores.tmp, na.rm = TRUE), 
                         3)