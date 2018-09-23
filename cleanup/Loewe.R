Loewe = function (response.mat,correction = T, Emin = 0, Emax = NA,nan.handle = c("L4")) {
  # scores <- list()
  #  method <- "Loewe"
  out <- tryCatch(Loewe2(response.mat,correction = correction, Emin = Emin, Emax = Emax, 
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

Loewe2 = function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
                   nan.handle = c("L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD2(response.mat, Emin = Emin, 
                                          Emax = Emax, nan.handle = nan.handle)$corrected.mat
  }
  single.fit <- FittingSingleDrug2(response.mat,fixed = c(NA, NA, NA, NA), nan.handle = nan.handle)
  #  single.fit <- 
  drug.col.model <- single.fit$drug.col.model
  drug.col.par <- coef(drug.col.model)
  d1.fun <- function(conc, drug.col.model) {
    drug.col.par <- coef(drug.col.model)
    # LL.4
    if(length(grep("LL.4", drug.col.model$call$fct))> 0 )
      (drug.col.par[3] + drug.col.par[2] * (conc/drug.col.par[4])^drug.col.par[1])/(1 + 
                                                                                      (conc/drug.col.par[4])^drug.col.par[1])
    else # L.4
      (drug.col.par[2] + (drug.col.par[3] - drug.col.par[2])/(1+exp(drug.col.par[1]*(conc-drug.col.par[4]))))
  }
  
  drug.row.model <- single.fit$drug.row.model
  drug.row.par <- coef(drug.row.model)
  d2.fun <- function(conc, drug.row.model) {
    drug.row.par <- coef(drug.row.model)
    if(length(grep("LL.4", drug.row.model$call$fct))> 0 ) # LL.4
      (drug.row.par[3] + drug.row.par[2] * (conc/drug.row.par[4])^drug.row.par[1])/(1 + 
                                                                                      (conc/drug.row.par[4])^drug.row.par[1])
    else #L.4
      (drug.row.par[2] + (drug.row.par[3] - drug.row.par[2])/(1+exp(drug.row.par[1]*(conc-drug.row.par[4]))))
  }
  row.conc <- as.numeric(rownames(response.mat))[-1]
  col.conc <- as.numeric(colnames(response.mat))[-1]
  loewe.mat <- response.mat
  
  # Four functions to calculate loewe 
  eq.LL4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {# Eq.8 in the ZIP paper
    x1/(drug.col.par[4] * (((x - drug.col.par[3])/(drug.col.par[2] - x))^(1/drug.col.par[1]))) + 
      x2/(drug.row.par[4] * (((x - drug.row.par[3])/(drug.row.par[2] - x))^(1/drug.row.par[1]))) - 1
  }
  
  eq.L4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
    x1/(drug.col.par[4] + log((drug.col.par[3]-x)/(x-drug.col.par[2]))/drug.col.par[1]) +
      x2/(drug.row.par[4] + log((drug.row.par[3]-x)/(x-drug.row.par[2]))/drug.row.par[1]) -1
  }
  
  eq.LL4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
    x1/(drug.col.par[4] * (((x - drug.col.par[3])/(drug.col.par[2] - x))^(1/drug.col.par[1]))) +
      x2/(drug.row.par[4] + log((drug.row.par[3]-x)/(x-drug.row.par[2]))/drug.row.par[1]) -1
  }
  
  eq.L4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
    x1/(drug.col.par[4] + log((drug.col.par[3]-x)/(x-drug.col.par[2]))/drug.col.par[1]) +
      x2/(drug.row.par[4] * (((x - drug.row.par[3])/(drug.row.par[2] - x))^(1/drug.row.par[1]))) - 1
  }
  
  cond1 = length(grep("LL.4", drug.col.model$call$fct))> 0 
  cond2 = length(grep("LL.4", drug.row.model$call$fct))> 0 
  if( cond1 == T & cond2 == T) eq = eq.LL4.LL4
  if( cond1 == T & cond2 == F) eq = eq.LL4.L4
  if( cond1 == F & cond2 == T) eq = eq.L4.LL4
  if( cond1 == F & cond2 == F) eq = eq.L4.L4
  
  for (i in 1:length(col.conc)) {
    for (j in 1:length(row.conc)) {
      x1 <- col.conc[i]
      x2 <- row.conc[j]
      
      options(warn = -1)
      slv = tryCatch(
        {
          slv <- nleqslv(max(drug.col.par[2] + 1, drug.row.par[2] + 1), eq, method = "Newton", x1=x1, x2=x2, drug.col.par = drug.col.par, drug.row.par = drug.row.par)
          
        }, error = function(cond){
          slv = list(termcd = 0, x = 0)
        }
      )
      
      
      options(warn = 0)
      
      if (slv$termcd == 1) {
        loewe.mat[j + 1, i + 1] <- slv$x
      }
      else {
        y.loewe1 <- d1.fun(x1 + x2, drug.col.model)
        y.loewe2 <- d2.fun(x1 + x2, drug.row.model)
        loewe.mat[j + 1, i + 1] <- ifelse(100 > max(y.loewe1, y.loewe2), max(y.loewe1, y.loewe2), 100)
      }
    }
  }
  
  return(response.mat - loewe.mat)
}
