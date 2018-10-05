ZIP2 = function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
          nan.handle = c("L4")) 
{
  if (correction) {
    nan.handle <- match.arg(nan.handle)
    response.mat <- BaselineCorrectionSD2(response.mat, NA, NA, nan.handle)$corrected.mat
  }
  single.fitted <- FittingSingleDrug2(response.mat, fixed = c(NA, Emin, Emax, NA), nan.handle)
  drug.col.response <- single.fitted$drug.col.fitted
  drug.row.response <- single.fitted$drug.row.fitted
  updated.single.mat <- mat.or.vec(nrow(response.mat), ncol(response.mat))
  colnames(updated.single.mat) <- colnames(response.mat)
  rownames(updated.single.mat) <- rownames(response.mat)
  
  updated.single.mat[1, c(2:ncol(response.mat))] <- drug.col.response[-1]
  updated.single.mat[c(2:nrow(response.mat)), 1] <- drug.row.response[-1]
  updated.col.mat <- updated.single.mat
  
  for (i in 2:ncol(response.mat)) {
    tmp <- as.data.frame(mat.or.vec(nrow(response.mat) - 1, 0))
    tmp$dose <- as.numeric(rownames(response.mat))[-1]
    tmp$dose[which(tmp$dose==0)] =  10^-10 # nonzero concentrations to take the log
    tmp$inhibition <- response.mat[c(2:nrow(response.mat)), i]
    
    tmp.min <- updated.single.mat[1, i]
    if (var(tmp$inhibition, na.rm = TRUE) == 0 | length(tmp$inhibition) == 1) {
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }
    
    if(nrow(tmp) == 1){
      # # no fitting
      fitted.inhibition = response.mat[c(2:nrow(response.mat)), i]
    }else{
      tmp.model = tryCatch(
        {
          tmp.model <- drm(inhibition ~ dose, data = tmp, fct = LL.4(fixed = c(NA, tmp.min, Emax, NA)), na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
        }, 
        warning = function(w) {
         tmp.model <- drm(inhibition ~ log(dose), data = tmp, fct = L.4(fixed = c(NA, tmp.min, Emax, NA)),
                          na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
        },
        error = function(e) {
          tmp.model <- drm(inhibition ~ log(dose), data = tmp, fct = L.4(fixed = c(NA, tmp.min, Emax, NA)), 
                           na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
          
        }
      )
      fitted.inhibition = suppressWarnings(fitted(tmp.model))
    }
    
    tmp$fitted.inhibition = fitted.inhibition
    # if (tmp$fitted.inhibition[nrow(response.mat) - 1] < 0) 
    #  tmp$fitted.inhibition[nrow(response.mat) - 1] <- tmp.min
    updated.col.mat[c(2:nrow(response.mat)), i] <- tmp$fitted.inhibition
  }
  
  updated.row.mat <- updated.single.mat
  for (i in 2:nrow(response.mat)) {
    tmp <- as.data.frame(mat.or.vec(ncol(response.mat) - 1, 0))
    tmp$dose <- as.numeric(colnames(response.mat))[-1]
    tmp$dose[which(tmp$dose==0)] =  10^-10 # nonzero concentrations to take the log
    
    tmp$inhibition <- response.mat[i, c(2:ncol(response.mat))]
    tmp.min <- updated.single.mat[i, 1]
    
    if (var(tmp$inhibition, na.rm = TRUE) == 0 | length(tmp$inhibition) == 1) {
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }
    
    if(nrow(tmp)==1){
      fitted.inhibition = response.mat[i, c(2:ncol(response.mat))]
    }else{
      tmp.model = tryCatch(
        {
          tmp.model <- drm(inhibition ~ dose, data = tmp, fct = LL.4(fixed = c(NA, tmp.min, Emax, NA)), na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
        }, 
        warning = function(w) {
          tmp.model <- drm(inhibition ~ log(dose), data = tmp, fct = L.4(fixed = c(NA, tmp.min, Emax, NA)), 
                           na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
          
        },
        error = function(w) {
          tmp.model<-drm(inhibition ~ log(dose), data = tmp, fct = L.4(fixed = c(NA, tmp.min, Emax, NA)), 
                         na.action = na.omit, control = drmc(errorm = FALSE, noMessage = T))
        }
        
      )
      fitted.inhibition = suppressWarnings(fitted(tmp.model))
    }
    
    # if (tmp$fitted.inhibition[ncol(response.mat) - 1] < 0) 
    #  tmp$fitted.inhibition[ncol(response.mat) - 1] <- tmp.min
    tmp$fitted.inhibition = fitted.inhibition
    updated.row.mat[i, c(2:ncol(response.mat))] <- tmp$fitted.inhibition
  }
  
  fitted.mat <- (updated.col.mat + updated.row.mat)/2
  zip.mat <- updated.single.mat
  for (i in 2:nrow(updated.single.mat)) {
    for (j in 2:ncol((updated.single.mat))) {
      zip.mat[i, j] <- updated.single.mat[i, 1] + updated.single.mat[1, j] - updated.single.mat[i, 1] * updated.single.mat[1, j]/100
    }
  }
  
  # fitted.mat[1, 1] <- 0
  # zip.mat[1, 1] <- 0
  # fitted.mat <- apply(fitted.mat, c(1, 2), function(x) ifelse(x > 100, 100, x))
  
  delta.mat <- (fitted.mat - zip.mat)
  delta.mat
}
