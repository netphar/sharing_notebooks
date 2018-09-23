CalculateSynergy = function (data, method = "ZIP", correction = TRUE, Emin = 0, Emax = NA, nan.handle = c("L4"))
{
  if (!is.list(data)) {
    stop("Input data is not a list format!")
  }
  if (!method %in% c("ZIP", "HSA", "Bliss", "Loewe")) {
    stop("The method parameter can only be one of the following: ZIP, HSA, Bliss and Loewe.")
  }
  dose.response.mats <- data$dose.response.mats
  num.pairs <- length(dose.response.mats)
  scores <- list()
  nan.handle <- match.arg(nan.handle)
  pb <- txtProgressBar(min = 0, max = num.pairs, style = 3)
  
  for (i in 1:num.pairs) {
    setTxtProgressBar(pb, i)
    
    response.mat <- dose.response.mats[[i]]
    scores[[i]] <- switch(method, ZIP = ZIP(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle), 
                          HSA = HSA(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle), 
                          Bliss = Bliss(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle), 
                          Loewe = Loewe(response.mat, correction, Emin = Emin, Emax = Emax, nan.handle = nan.handle))
  }
  data$scores <- scores
  data$method <- method
  close(pb)
  
  return(data)
}
