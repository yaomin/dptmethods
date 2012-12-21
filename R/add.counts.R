add.counts <-
function(x, ...) {
  ## if(nargs() < 3) add.counts.2smp.v2(x, ...)
  ## else add.counts.2smp.v2(x, add.counts.2smp.v2(...))
  res <- join.samples(x, ...)
  res <- as.data.frame(cbind(res[,1], rowSums(res[,-1])))
  names(res) <- c('poi', 'c1')
  res
}
