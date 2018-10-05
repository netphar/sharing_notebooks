own_rank = function(x){
  x_unique <- unique(x)
  x_ranks <- rank(x_unique)
  x_ranks[match(x,x_unique)]}