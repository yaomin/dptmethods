event.prob3 <-
function(x, sample.label,
         event, p, fun="quantile",p.fun="median",probs=0.27)
  ## assume x is in w3.res format as above
  ## p is the p3
{
    out <- NULL
    if(length(probs)==1) probs=rep(probs, length(event))
    event.ext <- dptmethods:::pattern.design(event, sample.label)
    for (i in seq(along=event)) {
	if(event[i] != 1) {
	    x[,i] <- 1-x[,i]
	    p[i] <- 1-p[i]
	}
	.P <- t(apply(x, 
		      1, 
		      function(x,sample.label) tapply(x, 
						      sample.label, 
						      fun, 
						      probs=probs[i]),
		      sample.label=sample.label))
	.p <- tapply(p, sample.label, p.fun)

	out <- c(out, list(apply(t(apply(.P,1, function(x) x*.p)), 1, prod)))
    }
    out <- as.data.frame(out)
    names(out)=paste("P", apply(event, 2, paste, collapse=""),sep="")
    out.rs <- rowSums(out) 
    out/out.rs  
}
