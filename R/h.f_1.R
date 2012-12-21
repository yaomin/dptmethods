H.f <-
function(theta, p) {
  theta*log(theta/p)+(1-theta)*log((1-theta)/(1-p))
}
