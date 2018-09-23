correctingConc = function(x,y,z) {
  
  y <- as.numeric(y)
  z <- as.numeric(z)
  
  #we choose rows where there is a correct combo of NSC1&2 id's. NSC is an internal drug identifier, 
  #here is the translation table https://wiki.nci.nih.gov/download/attachments/338237347/ComboCompoundNames_small.txt?version=1&modificationDate=1493822467000&api=v2)
  #NA for NSC2 column means that it was a single dose test
  x <- dplyr::filter(x, ((NSC1 == y & NSC2 == z) | (NSC1 == y & is.na(NSC2)) | (NSC1 == z & is.na(NSC2)) ))
  
  x_comb <- list()
  x_comb <- dplyr::filter(x, ( (NSC1 == y & NSC2 == z) ) )
  # this groups by column values to calculate the mean percent growth
  x_comb <- x_comb %>%
    dplyr::group_by(NSC1, NSC2, CONC1, CONC2) %>%
      dplyr::mutate(mean_single = mean(PERCENTGROWTH)) %>% 
        dplyr::arrange()
  
  
  aCONC11 <- unique(dplyr::filter(x_comb, (!is.na(CONC2))  )$CONC1)
  aNSC11 <- y
  aCONC21 <- unique(dplyr::filter(x_comb, (!is.na(CONC2))  )$CONC2)
  aNSC21 <-z
  
  x_one <- filter(x, (CONC1 %in% aCONC11)& NSC1 ==aNSC11 & is.na(CONC2))
  x_one <- x_one %>%
    dplyr::group_by(CONC1) %>%
      dplyr::mutate(mean_single = mean(PERCENTGROWTH)) %>% 
        dplyr::arrange()
  df <- data.frame()
  df <- as.data.frame(cbind(order(aCONC11), aCONC11))
  x_one <- dplyr::inner_join(x_one, df, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC11"))
  x_one <- x_one[!duplicated(x_one[c("CONC1")]),]
  
  x_two <- dplyr::filter(x, (CONC1 %in% aCONC21)& NSC1 ==aNSC21 & is.na(CONC2))
  x_two <- x_two %>%
    dplyr::group_by(CONC1) %>%
      dplyr::mutate(mean_single = mean(PERCENTGROWTH)) %>% # this gets mean values for all the experiments where one drug is tested by itself in cell line == 786-0
        dplyr::arrange()
  df2 <- data.frame()
  df2 <- as.data.frame(cbind(order(aCONC21), aCONC21))
  x_two <- dplyr::inner_join(x_two, df2, by=c("CONCINDEX1" = "V1", "CONC1" = "aCONC21"))
  x_two <- x_two[!duplicated(x_two[c("CONC1")]),]
  
  x_comb <- rbind(x_comb, x_one)
  x_comb <- rbind(x_comb, x_two)
  return(x_comb)
}