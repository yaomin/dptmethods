est.NB <-
function(tb) {
  lambda <- est.pois.lambda(tb)
  with(tb, freq[count==1])/dpois(1, lambda)
}
