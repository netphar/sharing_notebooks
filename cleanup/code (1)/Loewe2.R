Loewe2 = function (response.mat, correction = TRUE, Emin = NA, Emax = NA, 
          nan.handle = c("L4")) 
{
  if (correction) {
    response.mat <- BaselineCorrectionSD2(response.mat, Emin = Emin, 
                                         Emax = Emax, nan.handle)$corrected.mat
  }
  single.fit <- FittingSingleDrug2(response.mat,fixed = c(NA, NA, NA, NA), nan.handle)
  drug.col.model <- single.fit$drug.col.model
  drug.col.par <- coef(drug.col.model)
  
  drug.row.model <- single.fit$drug.row.model
  drug.row.par <- coef(drug.row.model)
  
  
  row.conc <- as.numeric(rownames(response.mat))[-1]
  col.conc <- as.numeric(colnames(response.mat))[-1]
  
  row.conc[which(row.conc==0)] = 10^-10 # avoid log(0)
  col.conc[which(col.conc==0)] = 10^-10 # avoid log(0)
  loewe.mat <- response.mat
  
  cond1 = length(grep("LL.4", drug.col.model$call$fct))> 0 
  cond2 = length(grep("LL.4", drug.row.model$call$fct))> 0 
  if( cond1 == T & cond2 == T) {eq = eq.LL4.LL4}
  if( cond1 == T & cond2 == F) {eq = eq.LL4.L4}
  if( cond1 == F & cond2 == T) {eq = eq.L4.LL4}
  if( cond1 == F & cond2 == F) {eq = eq.L4.L4}
  
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
      
      if (slv$termcd == 1) {
        loewe.mat[j + 1, i + 1] <- slv$x
      }
      else {
        y.loewe1 <- d1.fun(x1, x2, drug.col.model) #x1 col, x2 row
        y.loewe2 <- d2.fun(x1, x2, drug.row.model) #x1 col, x2 row
        loewe.mat[j + 1, i + 1] <- ifelse(100 > max(y.loewe1, y.loewe2), max(y.loewe1, y.loewe2), 100)
      }
    }
  }

  return(response.mat - loewe.mat)
}
