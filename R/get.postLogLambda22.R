get.postLogLambda22 <-
function(x, log=F) {
  if(log) {
    res <- log(x$paraCRM$lambda+0)*x$w.win[,2]^2
  }
  else res <- x$paraCRM$lambda*x$w.win[,2]^2
}
