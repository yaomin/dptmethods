event.v.prob <-
function(x,v,  event=c(0,0,0,1,1,1,1,1,1)) {
  p <- x
  for (i in seq(event)) {
    if(event[i] !=1) p[,i] <-  1-x[,i]
  }
  logp <- log(p)
  logw <- v/p^2
  logres <- logp*logw
  res <- apply(logres,1,sum)
  res
}
