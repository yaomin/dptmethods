#' @export
###. Initial filtering
##____________________________________________________________________
est.pois.lambda <- function(tab) {
  f1 <- with(tab, freq[count==1])
  f2 <- with(tab, freq[count==2])
  2*f2/f1
}

est.fdr <- function(n,tb) {
  lambda <- est.pois.lambda(tb)
  n.total <- with(tb, freq[count==1])/dpois(1, lambda)
  n.total*dpois(n,lambda)/with(tb, freq[count == n])
}

compute.initfilter.cutoff <- function(tbs, cutoff=0.5) {
  ## which count cutoff should be used so that at least one sample start to contain at least 'cutoff'
  ##  positives (max of counts across all samples should be at east this number)
  ifsatisfy <- rowSums(sapply(tbs, function(x) sapply(1:10, est.fdr, tb=x)) < cutoff) >0
  which(ifsatisfy)[1]
}

get.initfilter.cutoff <- function(counts.tab, cutoff=0.5, bysample=FALSE) {
  ##get.mappedCountTab()
  if(!bysample) {
    compute.initfilter.cutoff(counts.tab, cutoff=cutoff)
  } else {
    # Returns a vector with a per-sample cutoff rather than just the maximum cutoff
    ifsatisfy <- sapply(counts.tab, function(x) sapply(1:10, est.fdr, tb=x)) < cutoff
    apply(ifsatisfy, MARGIN=2, FUN=function(x) which(x)[1])
  }  
}
