normalize.sig <-
function(sig, bg) {
  centr <- apply(bg, 2, function(x) mean(x))
  diff <- centr-mean(centr)
  scl <- apply(bg, 2, function(x) sd(x))
  rescl <- scl/mean(scl)
  ##browser()
  for(i in seq(ncol(sig))) {
    sig[,i] <- (sig[,i]+diff[i])/rescl[i]
  }
  sig
}
