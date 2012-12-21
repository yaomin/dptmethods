extend.tbs <-
function(tbs) {
  lapply(tbs, function(x) rbind(c(0,est.NB(x)*dpois(0, est.pois.lambda(x))), x))  
}
