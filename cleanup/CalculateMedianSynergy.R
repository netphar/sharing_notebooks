CalculateMedianSynergy <- function(x) {
  # this adds an extra nested list to the calculated synergy values list. It contains median synergy rounded to 3 sig digits, like in synergyfinder
  
  x$median <- list()
  for (i in 1:nrow(x$drug.pairs)) {
    ifelse(is.na(x$scores[[i]]),
           x$median[[i]] <- NA, 
           x$median[[i]] <- round(mean(x$scores[[i]]), digits = 5) )
  }
  return(x)
}