###. Initial filtering
##____________________________________________________________________
est.pois.lambda <- function(tab) {
  f1 <- with(tab, freq[count==1])
  f2 <- with(tab, freq[count==2])
  2*f2/f1
}

# est.fdr <- function(n,tb) {
#   lambda <- est.pois.lambda(tb)
#   n.total <- with(tb, freq[count==1])/dpois(1, lambda)
#   n.total*dpois(n,lambda)/with(tb, freq[count == n])
# }

est.fdr <- function(n, tb) {
  ##browser()
  lambda.0 <- est.pois.lambda(tb)
  n0 <-  with(tb, freq[count==1]/dpois(1, lambda.0))
  cumcnt <- sapply(n, function(x) with(tb, sum(freq[as.integer(count)>x]))) 
  (ppois(n, lambda.0, lower=F)*n0)/cumcnt
}

compute.initfilter.cutoff <- function(tbs, cutoff=0.10) {
  ## which count cutoff should be used so that at least one sample start to contain at least 'cutoff'
  ##  positives (max of counts across all samples should be at east this number)
  ifsatisfy <- rowSums(sapply(tbs, function(x) sapply(1:10, est.fdr, tb=x)) < cutoff) ==length(tbs)
  which(ifsatisfy)[1]
}

get.initfilter.cutoff <- function(counts.tab, cutoff=0.10) {
  ##get.mappedCountTab()
  compute.initfilter.cutoff(counts.tab, cutoff=cutoff)
}