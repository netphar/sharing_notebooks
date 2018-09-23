Bliss = function (response.mat,correction = T, Emin = 0, Emax = NA,nan.handle = c("L4")) {
  # scores <- list()
  #  method <- "ZIP"
  out <- tryCatch(Bliss2(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
                         nan.handle = nan.handle), error = function(e) NA)
  #  {
  #  x$method <- method
  #  x$scores <- x$dose.response.mats
  #  x$scores[!is.na(x$scores)] <- NA
  #  x$scores <- replace(x$scores, !is.na(x$scores), NA)
  #  purrr::map(!is.na(x$scores), NA)
  #  x$scores <- lapply(x$scores, function(z) ifelse(!is.na(z), NA, z))
  #  return(x)
  #  })
  return(out)
}

Bliss2 = function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
                   nan.handle = c("L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD2(response.mat, Emin = Emin, 
                                          Emax = Emax, nan.handle)$corrected.mat
  }
  drug1.response <- response.mat[, 1]
  drug2.response <- response.mat[1, ]
  ref.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      ref.mat[i, j] <- drug1.response[i] + drug2.response[j] - 
        drug1.response[i] * drug2.response[j]/100
    }
  }
  syn.mat <- response.mat - ref.mat
  syn.mat
}

