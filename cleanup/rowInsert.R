rowInsert = function(x, blockID){
  class_old <- list()
  class_old <- sapply(x, class)
  #  print(class_old)
  NSC2_un <- vector(mode="numeric", length=1)
  NSC2_un <- unique(x$NSC2[!is.na(x$NSC2)])
  
  x <- as.tibble(t(apply(x, 1, exchange, NSC2_un)))
  x[] <- Map(`class<-`,x, class_old) #brilliant! https://stackoverflow.com/questions/27435873/changing-class-of-data-frame-columns-using-strings
  
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
  x <- x[-match(c("COMBODRUGSEQ","SCREENER","STUDY", "TESTDATE", "PLATE", "PANELNBR", "CELLNBR", "PREFIX1", "NSC1", "SAMPLE1", "CONCINDEX1", "CONC1", "CONCUNIT1", "PREFIX2", "NSC2", "SAMPLE2", "CONCINDEX2", "CONC2", "CONCUNIT2", "PERCENTGROWTH", "PERCENTGROWTHNOTZ","TESTVALUE", "CONTROLVALUE", "TZVALUE", "EXPECTEDGROWTH", "SCORE", "VALID", "PANEL", "CELLNAME", "mean_single"), names(x))]
  return(x)
}