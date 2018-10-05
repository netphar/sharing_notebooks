d1.fun <- function(col_conc, row_conc, drug.col.model) {
  drug.col.par <- coef(drug.col.model)
  # LL.4, conc must be raw
  if(length(grep("LL.4", drug.col.model$call$fct))> 0 ){
    conc = col_conc + row_conc
    (drug.col.par[3] + drug.col.par[2] * (conc/drug.col.par[4])^drug.col.par[1])/(1 + (conc/drug.col.par[4])^drug.col.par[1])
  } else # L.4, conc must be logscaled, ie. log(conc)
  {
    conc = log(col_conc+row_conc)
    (drug.col.par[2] + (drug.col.par[3] - drug.col.par[2])/(1+exp(drug.col.par[1]*(conc-drug.col.par[4]))))      
  }
  
}