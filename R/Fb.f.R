Fb.f <-
function(k,m,p) {
  res <- rep(NA, length(k))
  res[k<0] <- 0
  res[k>=0] <- pbinom(k[k>=0], m, p)
  res
}
