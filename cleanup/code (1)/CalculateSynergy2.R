CalculateSynergy2 = function (data, method = "ZIP", correction = TRUE, Emin = NA, 
          Emax = NA, nan.handle = c("L4")) 
{
  if (!is.list(data)) {
    stop("Input data is not a list format!")
  }
  if (!method %in% c("ZIP", "HSA", "BLISS", "LOEWE")) {
    stop("The method parameter can only be one of the following: ZIP, HSA, Bliss and Loewe.")
  }
  dose.response.mats <- data$dose.response.mats
  num.pairs <- length(dose.response.mats)
  scores <- list()
  nan.handle <- match.arg(nan.handle)
  for (i in 1:num.pairs) {
    response.mat <- dose.response.mats[[i]]
    
    # impute missing values
    # missing value imputation
    missing_index = which(is.na(response.mat), arr.ind = T)
    if(length(missing_index) !=0 ){
      for(i in 1:nrow(missing_index)){
        r = missing_index[i,1]
        c = missing_index[i,2]
        
        tmp = mean(c(response.mat[r+1,c],response.mat[r-1,c],response.mat[r,c-1],response.mat[r,c+1]),na.rm = T)
        if (is.na(tmp)) tmp = 0 # if no neighbors are found
        response.mat[r,c] = tmp
        
      }
    }
    
    scores[[i]] <- switch(method, ZIP = ZIP2(response.mat, 
                                            correction, Emin = Emin, Emax = Emax, nan.handle), 
                          HSA = HSA2(response.mat, correction, Emin = Emin, 
                                    Emax = Emax, nan.handle), BLISS = Bliss2(response.mat, 
                                                                            correction, Emin = Emin, Emax = Emax, nan.handle), 
                          LOEWE = Loewe2(response.mat, correction, Emin = Emin, 
                                        Emax = Emax, nan.handle))
  }
  data$scores <- scores
  data$method <- method
  return(data)
}