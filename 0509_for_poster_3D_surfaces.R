#no housekeeping. Put it

library('synergyfinder')

setwd('/Users/zagidull/Documents/fimm_files/synergy_calc_august/from_jing_dropbox')

read.csv(file = 'response_oneil_final-2.csv') -> oneil_final

#to get max response oneil
oneil_final[which(oneil_final$response == '100.076'),]
oneil_final[which(oneil_final$cell_line_name == 'A2780' & oneil_final$drug_row == 'L778123' & oneil_final$drug_col == 'BORTEZOMIB'),]

#loading bliss from NCI
readRDS('/Users/zagidull/Documents/fimm_files/synergy_calc_august/0309_reshaped.Bliss') -> reshaped_bliss

#adds new list with median synergy
CalculateMedianSynergy <- function(x) {
  x$median <- list()
  for (i in 1:nrow(x$drug.pairs)) {
    ifelse(is.na(x$scores[[i]]),
           x$median[[i]] <- NA, 
           x$median[[i]] <- round(mean(x$scores[[i]]), digits = 3) )
  }
  return(x)
}

CalculateMedianSynergy(reshaped_bliss) -> reshaped_with_mean_synergy

#get top 5 max synergies from 
head(sort(unlist(lapply(reshaped_with_mean_synergy$median,FUN=max)), decreasing = T), 5)
#get min synergy

#which block corresponds to the selected synergy
#which(reshaped_with_mean_synergy$median == as.numeric(57.597)) -> mx
#which(reshaped_with_mean_synergy$median == as.numeric(-42.274)) -> mn
#
mx <- c(8,22,23)
graph <- list()
graph$dose.response.mats <- reshaped_with_mean_synergy$dose.response.mats[mx]
graph$drug.pairs <- reshaped_with_mean_synergy$drug.pairs[mx,]
graph$scores <- reshaped_with_mean_synergy$scores[mx]
graph$method <- reshaped_with_mean_synergy$method

graph$drug.pairs$drug.row <- c('Lomustine', 'Vinorelbine', 'Vinorelbine')
rep('Thioguanine', 3) -> graph$drug.pairs$drug.col

synergyfinder::PlotSynergy(graph,type = '3D', save.file = F )


8 79037      752 Lomustine Thioguanine ADR-RES

22   608210 Vinorelbine      752        uM        uM 14:NCI/ADR-RES (High grade ovarian serous adenocarcinoma)
23   608210      752        uM        uM      14:SW-620  (Colon adenocarcinoma, lymph node)

#for 3D surfaces
graph <- list()
graph$dose.response.mats <- reshaped_with_mean_synergy$dose.response.mats[mx]
graph$dose.response.mats[[2]] <- reshaped_with_mean_synergy$dose.response.mats[[mn]]
graph$drug.pairs <- reshaped_with_mean_synergy$drug.pairs[mx,]
graph$drug.pairs[2,] <- reshaped_with_mean_synergy$drug.pairs[mn,]
graph$scores[[1]] <- reshaped_with_mean_synergy$scores[[mx]]
graph$scores[[2]] <- reshaped_with_mean_synergy$scores[[mn]]
graph$method <- reshaped_with_mean_synergy$method

#we need correct drug names
# source https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-ALMANAC?preview=/338237347/347474319/ComboCompoundNames_all.txt
# 85998 = Streptozotocin
# 141540 = Etoposide
# 218321 = Pentostatin
# 721517 = Zoledronate
#graph$drug.pairs$drug.row <- c('Streptozotocin', 'Etoposide')
#graph$drug.pairs$drug.col <- c('Pentostatin', 'Zoledronate')

