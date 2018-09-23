CalculateSynergy2 = function (data, method = "ZIP", correction = TRUE, Emin = 0, 
          Emax = 100, nan.handle = c("L4")) 
{
  if (!is.list(data)) {
    stop("Input data is not a list format!")
  }
  if (!method %in% c("ZIP", "HSA", "Bliss", "LOEWE")) {
    stop("The method parameter can only be one of the following: ZIP, HSA, Bliss and Loewe.")
  }
  dose.response.mats <- data$dose.response.mats
  num.pairs <- length(dose.response.mats)
  scores <- list()
  nan.handle <- match.arg(nan.handle)
  for (i in 1:num.pairs) {
    response.mat <- dose.response.mats[[i]]
    scores[[i]] <- switch(method, ZIP = ZIP2(response.mat, 
                                            correction, Emin = Emin, Emax = Emax, nan.handle), 
                          HSA = HSA2(response.mat, correction, Emin = Emin, 
                                    Emax = Emax, nan.handle), Bliss = Bliss2(response.mat, 
                                                                            correction, Emin = Emin, Emax = Emax, nan.handle), 
                          LOEWE = Loewe2(response.mat, correction, Emin = Emin, 
                                        Emax = Emax, nan.handle))
  }
  data$scores <- scores
  data$method <- method
  return(data)
}