library('tidyverse')
library('synergyfinder')
library('reshape2')

#removing lists with na
A498[!is.na(A498)] -> A498_noNA

# bind all together
big_data_a498 = do.call(rbind, A498_noNA) 

#removing unneeded columns
big_data_a498 <- big_data_a498[-match(c("CONCUNIT2", "PERCENTGROWTH", "PERCENTGROWTHNOTZ","VALID","PANEL", "CELLNAME","mean_single"), names(big_data_a498))]

#counting columns and removing where it is not equal to 16. But care when we have 
big_data_a498 <- big_data_a498 %>%
  group_by(BlockID) %>%
    mutate(n = n()) %>%
      filter(n == 16) 

big_data_a498 <- big_data_a498[-match(c("n"), names(big_data_a498))]
big_data_test <- big_data_a498 %>%
  filter(DrugRow == 747972  & DrugCol == 715055)

715055
747972

fiver[-c(26:50),] -> fiver
try(ReshapeData(big_data_a498, data.type = 'inhibition'))

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