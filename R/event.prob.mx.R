event.prob.mx <-
function(x, event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9), k=1e-36)
  ## assume x is in w3.res format as above
  ## p is the p3
{
    
  res <- x
  for (i in seq(along=event)) {
    if(event[i] !=1) {
      x[,i] <-  1-x[,i]
      p[i] <- 1-p[i]
    }
  }

  for (i in seq(along=event)) {
    res[,i] <- log2((x[,i]+k)/(1-x[,i]+k)*((1-p[i])/p[i])) ## test on p + k
  }
  res
}
