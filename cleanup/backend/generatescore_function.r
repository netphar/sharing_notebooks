# Function to generate the response table with synergy scores

# response input must contain the following columns
# block_id    conc_r    conc_c    response    drug_row    drug_col    conc_r_unit    conc_c_unit    cell_line_name

# response = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\code\\response_template.csv", stringsAsFactors = F)
# setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\code")

GenerateScore = function(response){
  # add row and column numbers
  response = ddply(response, c("cell_line_name","drug_row","drug_col","block_id"), transform, row = own_rank(conc_r), col = own_rank(conc_c))
  scores = list()
  m = unique(response$block_id)
  
  # options(show.error.messages = F)
  for(i in m){
    cat('\r',i)
    index = which(response$block_id == i)
    data.tmp = response[index,]
    data.tmp2 = ReshapeData2(data.tmp, data.type = "inhibition")
    
    hsa.tmp = CalculateSynergy2(data.tmp2, method = "HSA", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
    bliss.tmp = CalculateSynergy2(data.tmp2, method = "BLISS", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
    zip.tmp = CalculateSynergy2(data.tmp2, method = "ZIP", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
    loewe.tmp = CalculateSynergy2(data.tmp2, method = "LOEWE", correction = T, nan.handle = "L4", Emin = NA, Emax = NA)
    
    data.tmp$synergy_zip = apply(data.tmp[,c("row","col")], 1, function(x) zip.tmp$scores[[1]][x[1],x[2]])
    data.tmp$synergy_hsa = apply(data.tmp[,c("row","col")], 1, function(x) hsa.tmp$scores[[1]][x[1],x[2]])
    data.tmp$synergy_bliss = apply(data.tmp[,c("row","col")], 1, function(x) bliss.tmp$scores[[1]][x[1],x[2]])
    data.tmp$synergy_loewe = apply(data.tmp[,c("row","col")], 1, function(x) loewe.tmp$scores[[1]][x[1],x[2]])
    
    scores[[i]] = data.tmp
    # if (unique(data.tmp$block_id)!=i) 
    #   print(i)
    # flush.console()
  }
  # options(show.error.messages = TRUE)
  
  response_with_scores = do.call(rbind, scores)
  write.csv(response_with_scores, "response_with_scores.csv", row.names = F)
  save.image("response_with_scores.RData")
  return(response_with_scores)

}


