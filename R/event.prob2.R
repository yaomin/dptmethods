event.prob2 <-
function(x, sample.label,
         event, p, fun="mean",...)
  ## assume x is in w3.res format as above
  ## p is the p3
{
  event.P <- t(apply(x, 
		     1, 
		     function(x,sample.label) tapply(x, 
						     sample.label, 
						     fun, 
						     ...),
		     sample.label=sample.label))
  event.p <- tapply(p, sample.label,  fun, ...)
 
  out <- NULL
  for (i in seq(ncol(event))) {
      .P <- event.P
      .p <- event.p
      for (j in seq(along=event[,i])) {
	  if(event[j,i] !=1) {
	      .P[,j] <- 1-.P[,j]
	      .p[j] <- 1-.p[j]
	  }
      }
      out <- c(out, list(apply(t(apply(.P,1, function(x) x*.p)), 1, prod)))
  }
  out <- as.data.frame(out)
  names(out)=paste("P", apply(event, 2, paste, collapse=""),sep="")
  out.rs <- rowSums(out) 
  out/out.rs  
}
