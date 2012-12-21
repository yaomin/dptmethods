compute.p <-
function(x, event=c(0,0,0,1,1,1,1,1,1))
  ## assume x is in w3.res format as above
{
  res <- x
  for (i in seq(event)) {
    if(event[i] !=1) res[,i] <-  1-x[,i]
  }
  res
}
