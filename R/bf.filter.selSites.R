bf.filter.selSites <-
function(sites.bf, cutoff,
                               which.score=c("score.1","score","score.2","score.P","qvalue")) {
  ##browser()
  which.score <- match.arg(which.score) 
  if(which.score %in% c("score.1","score","score.2")) {
    sites <- subset(sites.bf, subset=sites.bf[,which.score] > cutoff)
  } else if(which.score == "score.P") {
    sites <- subset(sites.bf, subset=sites.bf[,which.score] < cutoff)
  } else if(which.score == "qvalue") {
    pvalues <- sites.bf[,"score.P"]
    sites.qvalue <- rep(NA, length=length(pvalues))
    p.finite <- is.finite(pvalues)
    sites.qvalue[p.finite] <- compute.qvalues(pvalues[p.finite],lambda=0)
    sites <- data.frame(sites.bf, qvalue=sites.qvalue)
    sites <- subset(sites, subset= qvalue< cutoff)
  } else stop("Invalid bf score!")
  sites
}
