predicate = function(x,NSC2_un) {
  if (as.numeric(x[["CONCINDEX2"]]== 0) & as.numeric(x[["NSC1"]]) == NSC2_un) {
    return(T)
  }
  else {return(F)}
}