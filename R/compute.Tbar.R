compute.Tbar <-
function(T,v) {
  res <- T
  for (i in seq(ncol(T))) {
    res[,i] <- T[,i]*v[,i]
  }
  rowSums(res)/rowSums(v)
}
