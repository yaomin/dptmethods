denoiseSignal <-
function(tb) {
  ##browser()
  N <- sapply(tb, est.NB)
  lams <- sapply(tb, est.pois.lambda)
  newtb <- vector("list", length(tb))
  for (i in seq(along=tb)) {
    newtb[[i]] <- tb[[i]]
    newtb[[i]]$freq <- as.integer(tb[[i]]$freq-N[i]*dpois(tb[[i]]$count, lams[i]))
  }
  newtb
}
