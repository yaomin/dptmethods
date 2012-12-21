mean.qint <-
function(x, probs=c(1,0.75)) {
  q12 <- quantile(x, probs)
  weights <- x<=q12[1] & x>q12[2]
  weighted.mean(x, weights)
}
