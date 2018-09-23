# use the example data score.RData as the input;
smoothing = function (scores.dose, len = 2) 
{
  scores.dose = as.data.frame(scores.dose)
  options(scipen = 999)
  nr <- nrow(scores.dose)
  nc <- ncol(scores.dose)
  
  # missing value imputation
  missing_index = which(is.na(scores.dose),arr.ind = T)
  if(length(missing_index) !=0 ){
    for(i in 1:nrow(missing_index)){
      r = missing_index[i,1]
      c = missing_index[i,2]
      
      scores.dose[r,c] = mean(c(scores.dose[r+1,c],scores.dose[r-1,c],scores.dose[r,c-1],scores.dose[r,c+1]),na.rm = T) 
    }
  }
  
  extended.row.idx <- seq(1, nr, length = (nr - 1) * (len + 
                                                        2) - (nr - 2))
  extended.col.idx <- seq(1, nc, length = (nc - 1) * (len + 
                                                        2) - (nc - 2))
  krig.coord <- cbind(rep(extended.row.idx, each = length(extended.col.idx)), 
                      rep(extended.col.idx, length(extended.row.idx)))
  extended.scores <- SpatialExtremes::kriging(data = c(as.matrix(scores.dose)), 
                                              data.coord = cbind(rep(1:nr, nc), rep(1:nc, each = nr)), 
                                              krig.coord = krig.coord, cov.mod = "whitmat", grid = FALSE, 
                                              sill = 1, range = 10, smooth = 0.8)$krig.est
  extended.scores <- matrix(extended.scores, nrow = length(extended.row.idx), 
                            ncol = length(extended.col.idx), byrow = TRUE)
  
  # extended.scores = data.frame(extended.scores)
  extended.scores = round(extended.scores, 3)
  col.dose = as.numeric(colnames(scores.dose))
  row.dose = as.numeric(rownames(scores.dose))
  
  
  row.index = 1:nr
  col.index = 1:nc
  
  row.tick = unlist(lapply(row.index, function(x) which(extended.row.idx == x)))
  col.tick = unlist(lapply(col.index, function(x) which(extended.col.idx == x)))
  
  rownames(extended.scores)[row.tick] = row.dose
  colnames(extended.scores)[col.tick] = col.dose
  
  list1 = lapply(diff(row.dose)/diff(row.tick), function(x) cumsum(rep(x, each = len)))
  list2 = as.list(row.dose[-nr])
  rownames(extended.scores)[-row.tick] = as.character(round(unlist(mapply("+", list1, list2, SIMPLIFY = FALSE)),4))
  
  list1 = lapply(diff(col.dose)/diff(col.tick), function(x) cumsum(rep(x, each = len)))
  list2 = as.list(col.dose[-nc])
  colnames(extended.scores)[-col.tick] = as.character(round(unlist(mapply("+", list1, list2, SIMPLIFY = FALSE)),4))
  
  return(extended.scores)
  # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
  # heatmap.2(extended.scores, Rowv = F, Colv = F, dendrogram = 'none',trace = 'none',col = my_palette, density.info = "none")
  
}
