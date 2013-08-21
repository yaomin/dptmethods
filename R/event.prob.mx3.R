event.prob.mx3 <-
function(x, event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9))
  ## assume x is in w3.res format as above
  ## p is the p3
{
    
  for (i in seq(along=event)) {
    if(event[i] !=1) {
      x[,i] <-  1-x[,i]
      p[i] <- 1-p[i]
    }
  }
  list(x=x, p=p)
}
