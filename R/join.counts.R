join.counts <-
function(x, ...,by=100) {
  
  ends <- range.pois(x, ..., which.col=1)
  x.ext <- extend.seq(x, start=ends[1], end = ends[2],by=by)
  res <- lapply(list(...), function(x, ends) extend.seq(x,
                                                        start=ends[1],
                                                        end=ends[2],
                                                        by=by),
                ends = ends)
  res.cnt <- as.data.frame(sapply(res, function(x) x[,2]))
  res.cmb <- cbind(x.ext, res.cnt)
  names(res.cmb) <- c('poi', paste('c', seq(ncol(res.cmb)-1), sep=''))
  res.cmb
  
}
