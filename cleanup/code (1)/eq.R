# Four functions to calculate loewe 
eq.LL4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {# Eq.8 in the ZIP paper
  x1/(drug.col.par[4] * (((x - drug.col.par[3])/(drug.col.par[2] - x))^(1/drug.col.par[1]))) + 
    x2/(drug.row.par[4] * (((x - drug.row.par[3])/(drug.row.par[2] - x))^(1/drug.row.par[1]))) - 1
}

eq.L4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) { # x1, x2 to be log scaled
  log(x1)/(drug.col.par[4] + log((drug.col.par[3]-x)/(x-drug.col.par[2]))/drug.col.par[1]) +
    log(x2)/(drug.row.par[4] + log((drug.row.par[3]-x)/(x-drug.row.par[2]))/drug.row.par[1]) -1
}

eq.LL4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) { # x2 to be log-scaled
  x1/(drug.col.par[4] * (((x - drug.col.par[3])/(drug.col.par[2] - x))^(1/drug.col.par[1]))) +
    log(x2)/(drug.row.par[4] + log((drug.row.par[3]-x)/(x-drug.row.par[2]))/drug.row.par[1]) -1
}

eq.L4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) { # x1 to be log-scaled
  log(x1)/(drug.col.par[4] + log((drug.col.par[3]-x)/(x-drug.col.par[2]))/drug.col.par[1]) +
    x2/(drug.row.par[4] * (((x - drug.row.par[3])/(drug.row.par[2] - x))^(1/drug.row.par[1]))) - 1
}
