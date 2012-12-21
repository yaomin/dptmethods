join.samples <-
function(x,...,by=100) {
  
  ends <- range.pois(x, ..., which.col=2)
  x.ext <- extend.seq(x[,c(2,4)], start=ends[1], end=ends[2],by=by)
  if(nargs()>2) {
    res <- lapply(list(...),function(x, ends) extend.seq(x[,c(2,4)],
                                                         start=ends[1],
                                                         end=ends[2],
                                                         by=by),
                  ends=ends)
    res.cnt <- as.data.frame(sapply(res, function(x) x[,2]))
    res.cmb <- cbind(x.ext, res.cnt)
  } else {
    res.cmb <- x.ext
  }
  names(res.cmb) <- c('poi', paste('s', seq(ncol(res.cmb)-1), sep=''))
  res.cmb
  
}
