logit.f <-
function(x, k=1e-36) {
  log((x+k)/(1-x+k))
}
