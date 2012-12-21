consolidate.lows <-
function(dat, cols=2:3, limit=0, method=c("sum","max", "low")) {
  method <- match.arg(method)
  if(method == "sum") {
    idx <- rowSums(dat[,cols]) > limit
    
  }
  else if (method == "max") {
    idx <- apply(dat[, cols], 1, max) > limit
  }
  else if (method == "low") {
     idx <- apply(dat[, cols], 1, max) <= limit
  }
  else stop("Not implemented yet")
  res <- dat[idx,]
}
