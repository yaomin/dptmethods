R.merger <-
function(smp1,smp2) {
  if(missing(smp2)) {
    out <- smp1
    dimnames(out)[[2]] <- c("poi", "s1")
  } else {
    poi <- union(smp1[,1], smp2[,1])
    x1 <- matrix(0, nrow=length(poi), ncol=ncol(smp1)-1)
    x2 <- rep(0, length(poi))
    x1[match(smp1[,1], poi),] <- smp1[,-1]
    x2[match(smp2[,1], poi)] <- smp2[,2]
    out <- cbind(poi, x1, x2)
    dimnames(out)[[2]] <- c("poi", paste("s", 1:(ncol(x1)+1), sep=""))
  }
  out
}
