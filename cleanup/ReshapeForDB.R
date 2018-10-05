ReshapeForDB = function(x, y) {  
  # x is output of CalculateSynergy file. It is a list of lists, with dose.response.mats, scores, method, drug.pairs
  # y is the name of the Synergy method used for calculating the input file. In format SynergyHSA/ZIP/Bliss/Loewe
  
  # init number of iterations needed, find which blockID's have NA as synergy score and create a nested list of length iter to store results
  iter <- c(1:nrow(x$drug.pairs))
  to.populate <- which(is.na(x$scores))
  calculated.synergy.datalist <- vector('list', length(iter))
  
  temp.pb <- txtProgressBar(min = 0, max = length(iter), style = 3)
  
  for (i in iter) {
    
    # progress bar
    setTxtProgressBar(temp.pb, i)
    
    # precalc values of the drug.pairs. Useful if we need to add drug.row and drug.col. Otherwise probably unnecessary
    temp.info <- x$drug.pairs[i,]
    
    # get dose-response values in long matrix format
    temp.dose <- reshape2::melt(x$dose.response.mats[[i]])
    
    # get synergy scores in long
    # if score is NA, reconstruct a matrix with the same dimension as dose.response.mats filled with NA's
    if (i %in% to.populate) {
      temp.scores <- reshape2::melt(matrix(, nrow = nrow(x$dose.response.mats[[i]]), 
                                           ncol = ncol(x$dose.response.mats[[i]]), 
                                           dimnames = dimnames(x$dose.response.mats[[i]])))
    } else {
      temp.scores <- reshape2::melt(x$scores[[i]])
    }
    
    # merge both df's together
    temp <- base::merge(temp.dose,temp.scores, by = c('Var1', 'Var2'))
    
    # computer number of rows in long format of data in order to know how many times we need to paste the blockID
    nrows <- nrow(temp)
    
    # adjust names and add drug row/col and blockID's
    names(temp) <- c('ConcR','ConcC','ResponseInhibition', y)
    temp$BlockID <- rep(temp.info$blockIDs, times = nrows)
    temp$DrugRow <- rep(temp.info$drug.row, times = nrows)
    temp$DrugCol <- rep(temp.info$drug.col, times = nrows)
    temp$DrugColUnit <- rep(temp.info$concCUnit, times = nrows)
    temp$DrugRowUnit <- rep(temp.info$concRUnit, times = nrows)   
    calculated.synergy.datalist[[i]] <- temp
    
    # put smoothing.R here?
  }
  close(temp.pb)
  return(calculated.synergy.datalist)
  }
