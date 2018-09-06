readRDS('3008_big_data_with_Loewe') -> big_data_with_Loewe
readRDS('3008_datalist_with_Loewe') -> datalist_with_Loewe

cor(big_data_with_Loewe$Synergy_ZIP, big_data_with_Loewe$Synergy_Loewe, use = 'complete.obs')
cor(big_data_with_Loewe$Synergy_ZIP, big_data_with_Loewe$Synergy_Bliss, use = 'complete.obs')
cor(big_data_with_Loewe$Synergy_ZIP, big_data_with_Loewe$Synergy_HSA, use = 'complete.obs')
cor(big_data_with_Loewe$Synergy_Loewe, big_data_with_Loewe$Synergy_Bliss, use = 'complete.obs')
cor(big_data_with_Loewe$Synergy_Loewe, big_data_with_Loewe$Synergy_HSA, use = 'complete.obs')
cor(big_data_with_Loewe$Synergy_Bliss, big_data_with_Loewe$Synergy_HSA, use = 'complete.obs')

#remove zeroscor
cor(big_data_with_Loewe$Synergy_ZIP[which(big_data_with_Loewe$Synergy_ZIP != 0)], 
    big_data_with_Loewe$Synergy_Loewe[which(big_data_with_Loewe$Synergy_ZIP != 0)], use = 'complete.obs')

# pearson is about linear relationships and is calculated on true values
# spearman is about monotonicity and is calcualted on ranks 
# https://stats.stackexchange.com/a/14963
cor(big_data_with_Loewe$Synergy_ZIP, big_data_with_Loewe$Synergy_Loewe, use = 'complete.obs', method = 'spearman')

#spearman
#cor(big_data_with_Loewe$Synergy_ZIP, big_data_with_Loewe$Synergy_Loewe, use = 'complete.obs', method = 'spearman')
cor(big_data_with_Loewe$Synergy_ZIP, big_data_with_Loewe$Synergy_Bliss, use = 'complete.obs',method = 'spearman')
cor(big_data_with_Loewe$Synergy_ZIP, big_data_with_Loewe$Synergy_HSA, use = 'complete.obs',method = 'spearman')
cor(big_data_with_Loewe$Synergy_Loewe, big_data_with_Loewe$Synergy_Bliss, use = 'complete.obs',method = 'spearman')
cor(big_data_with_Loewe$Synergy_Loewe, big_data_with_Loewe$Synergy_HSA, use = 'complete.obs',method = 'spearman')
cor(big_data_with_Loewe$Synergy_Bliss, big_data_with_Loewe$Synergy_HSA, use = 'complete.obs',method = 'spearman')