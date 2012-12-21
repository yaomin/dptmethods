get.postMean.cutoffs.aux2 <-
function(pmeans, lambda, cut.p=0.05) {
  pm1 <- pmeans[pmeans<lambda]
  pm2 <- pmeans[pmeans>lambda]

  cut1 <- quantile(pm1, 1-cut.p)
  cut2 <- quantile(pm2, cut.p)
  (cut1+cut2)/2
}
