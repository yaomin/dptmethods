get.postLogLambda3 <-
function(x, log=F) {
  if(log) {
    res <- log(x$signal.win+0)*x$w.win[,3]
  }
  else res <- x$signal.win*x$w.win[,3]
}
