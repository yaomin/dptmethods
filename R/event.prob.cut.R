event.prob.cut <-
function(n, event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9), p.cut=0.001,k=1e-36) {
  res <- NULL
  for (i in seq(event)) {
    if(event[i] !=1) {
      p[i] <- 1-p[i]
    }
    res <- cbind(res, rbeta(n, p[i], 1-p[i]))
  }
  
  

  for (i in seq(event)) {
    res[,i] <- log2((res[,i]+k)/(1-res[,i]+k)*((1-p[i])/p[i])) ## test on p + k
  }
  quantile(apply(res, 1, min), 1-p.cut)
}
