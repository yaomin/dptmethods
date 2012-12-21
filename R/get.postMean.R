get.postMean <-
function(x, log=F) {
  if(log) {
    res <- log(x$signal.win+0)*x$w.win[,3]+log(x$paraCRM$lambda+0)*x$w.win[,2]
  }
  else res <- x$signal.win*x$w.win[,3]+x$paraCRM$lambda*x$w.win[,2]
  res
}
