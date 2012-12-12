## Filename: analysis.r
## By: Yaomin Xu (yaomin@gmail.com)
## Date: Tue Dec 11 16:28:38 2012
## Notes:   
## =============================================================================

###. count.tbs
##____________________________________________________________________
## Prostate cell line project folder
####  /proj/GenomicsHD/Ting_ProstateCell/analyses/Updated/DPT_ws/Preprocess
load("countdistr.RData")

cnt <- 0:15
fdrs <- sapply(count.tbs, function(x) sapply(cnt, est.fdr, tb=x))
(lambda.1 <- est.pois.lambda(count.tbs[[1]]))
(lambda.2 <- est.pois.lambda(count.tbs[[2]]))
(lambda.3 <- est.pois.lambda(count.tbs[[3]]))
## [1] 0.6540025
## [1] 0.5604019
## [1] 0.5909182

(N0.1 <- with(count.tbs[[1]], freq[count==1]/dpois(1, lambda.1)))
(N0.2 <- with(count.tbs[[2]], freq[count==1]/dpois(1, lambda.2)))
(N0.3 <- with(count.tbs[[3]], freq[count==1]/dpois(1, lambda.3)))


cal.count.tbs.0 <- function(tb) {
  ##browser()
  lambda.0 <- est.pois.lambda(tb)
  n0 <-  with(tb, freq[count==1]/dpois(1, lambda.0))
  freq0 <- as.integer(n0*dpois(as.numeric(tb$count), lambda.0))
  data.frame(count=tb$count, freq=freq0)
}

est.count <- function(k, tb) {
  sapply(k, function(x) with(tb, sum(freq[as.integer(count)>x])))
}

est.fdr2 <- function(n, tb) {
  ##browser()
  lambda.0 <- est.pois.lambda(tb)
  n0 <-  with(tb, freq[count==1]/dpois(1, lambda.0))
  cumcnt <- sapply(n, function(x) with(tb, sum(freq[as.integer(count)>x]))) 
  (ppois(n, lambda.0, lower=F)*n0)/cumcnt
}

count.tbs.0 <- list()
for(i in 1:3) {
  count.tbs.0[[i]] <- cal.count.tbs.0(count.tbs[[i]])
}

prop.tbs.0 <- list()
for(i in 1:3) {
  prop0 <- count.tbs.0[[i]][,2]/count.tbs[[i]][,2]
  p0 <- 
  prop.tbs.0[[i]] <- data.frame(count=count.tbs.0[[i]]$count, prop=prop0)
  
}



pdf(file="initfilter.pdf", width=11, height=8.5)
with(count.tbs[[1]], plot(as.numeric(count),freq, type="l",col="red"))
with(count.tbs[[2]], lines(as.numeric(count),freq, type="l",col="blue"))
with(count.tbs[[3]], lines(as.numeric(count),freq, type="l",col="green"))


dev.off()
