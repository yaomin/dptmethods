bg.N <-
function(tbs) {
  bg.lambda <- mean(sapply(tbs, est.pois.lambda), trim=.1)
  bg.N <- median(sapply(tbs, est.NB))
  max.N <- max(sapply(tbs, function(x) max(x$count)))
  count <- seq(0, max.N)
  freq <- bg.N*dpois(count, bg.lambda)
  data.frame(count, freq)
}
