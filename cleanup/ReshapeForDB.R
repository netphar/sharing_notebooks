ReshapeForDB = function(x) {  
  # init number of iterations needed, find which blockID's have NA as synergy score and create a nested list of length iter to store results
  iter <- c(1:nrow(x$drug.pairs))
  to.populate <- which(is.na(x$scores))
  calculated.synergy.datalist <- vector('list', length(iter))
  
  for (i in iter) {
    
    # precalc values of the drug.pairs. Useful if we need to add drug.row and drug.col. Otherwise probably unnecessary
    temp.info <- x$drug.pairs[i,]
    
    # get dose-response values in long
    temp.dose <- melt(x$dose.response.mats[[i]])
    
    # get synergy scores in long
    # if score is NA, reconstruct a matrix with the same dimension as dose.response.mats filled with NA's
    if (i %in% to.populate) {
      temp.scores <- melt(matrix(, nrow = nrow(x$dose.response.mats[[i]]), ncol = ncol(x$dose.response.mats[[i]]), 
                                 dimnames = dimnames(x$dose.response.mats[[i]])))
    } else {
      temp.scores <- melt(x$scores[[i]])
    }
    
    # merge both df's together
    temp <- merge(temp.dose,temp.scores, by = c('Var1', 'Var2'))
    
    # computer number of rows in long format of data in order to know how many times we need to paste the blockID
    times <- nrow(temp)
    
    # adjust names and add drug row/col and blockID's
    names(temp) <- c('ConcR','ConcC','ResponseInhibition', 'SynergyHSA')
    temp$BlockID <- rep(temp.info$blockIDs, times = times)
    #    temp$DrugRow <- rep(temp.info$drug.row, times = times)
    #    temp$DrugCol <- rep(temp.info$drug.col, times = times)
    calculated.synergy.datalist[[i]] <- temp
    
  }
  return(calculated.synergy.datalist)
  }
