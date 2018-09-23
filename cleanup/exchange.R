exchange = function(x,NSC2_un) {
  predicate = function(x,NSC2_un) {
    if (as.numeric(x[["CONCINDEX2"]]== 0) & as.numeric(x[["NSC1"]]) == NSC2_un) {
      return(T)
    }
    else {return(F)}
  }
  if (predicate(x,NSC2_un)){
    placeholderIndex <- x[["CONCINDEX1"]]
    placeholderConc <- x[["CONC1"]]
    x[["CONCINDEX1"]] <- x[["CONCINDEX2"]]
    x[["CONCINDEX2"]] <- placeholderIndex
    x[["CONC1"]] <- x[["CONC2"]]
    x[["CONC2"]] <- placeholderConc
    return(x)
  } else { return(x)}
} 