#housekeeping
rm(list=ls(all=TRUE))

library('tidyverse')

setwd("/Users/zagidull/Documents/git/sharing_notebooks/cleanup")

Bliss.rbound <- readRDS('0110_rbound_Bliss')
HSA.rbound <- readRDS('0110_rbound_HSA')
Loewe.rbound <- readRDS('0110_rbound_Loewe')
ZIP.rbound <- readRDS('0110_rbound_ZIP')

Bliss.rbound %>% group_by(BlockID) %>% 

#one <- merge(Bliss.rbound, HSA.rbound, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID'))
#two <- merge(Loewe.rbound, ZIP.rbound, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID'))
#final <- merge(one, two, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID'))
  

#saveRDS(object = final, file = '0210_rbound_All')
  
  Bliss.rbound[Bliss.rbound$BlockID == '1:NCI/ADR-RES',] -> bliss
HSA.rbound[HSA.rbound$BlockID == '1:NCI/ADR-RES',] -> hsa
ZIP.rbound[ZIP.rbound$BlockID == '1:NCI/ADR-RES',] -> zip
Loewe.rbound[Loewe.rbound$BlockID == '1:NCI/ADR-RES',] -> loewe
  
out <-  merge(merge(ZIP.rbound, Loewe.rbound, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID')), 
        merge(HSA.rbound, Bliss.rbound, by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID')), 
        by = c('ConcR', 'ConcC', 'ResponseInhibition', 'BlockID'))

saveRDS(object = out, file = '0210_finishprocessing.R')
out1 <- out

out <- test


# input is original data, unsorted
input.all <- out1

out1 <- out

# arrange first by BlocKID, then ConcC, then ConcR
out1 <- out %>%  arrange(BlockID, ConcC, ConcR)

# to perform cor analysis
out <- out1

cor(out$SynergyZIP, out$SynergyLoewe, use = 'complete.obs') # [1] 0.4803294

cor(out$SynergyZIP, out$SynergyBliss, use = 'complete.obs') # [1] 0.9470076

cor(out$SynergyZIP, out$SynergyHSA, use = 'complete.obs') # [1] 0.9205561

cor(out$SynergyLoewe, out$SynergyBliss, use = 'complete.obs') # [1] 0.4698621

cor(out$SynergyLoewe, out$SynergyHSA, use = 'complete.obs') # [1] 0.4943656

cor(out$SynergyBliss, out$SynergyHSA, use = 'complete.obs') # [1] 0.9533624


out[out$BlockID=='1000:ACHN',] -> test


H B L Z

H:B H:L H:Z

B:L B:Z

L:Z

