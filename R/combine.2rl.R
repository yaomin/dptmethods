combine.2rl <-
function(x,y) {
  ## x and y are RangesList
  x.names <- names(x)
  y.names <- names(y)
  union.names <- union(x.names, y.names)
  com.names <- intersect(x.names,y.names)
  x.only.names <- setdiff(x.names, com.names)
  y.only.names <- setdiff(y.names, com.names)
  .outlist <- vector("list",length(union.names))
  names(.outlist) <- union.names
  if(length(com.names) >0) {
    for (nm in com.names) {
      .outlist[[nm]] <- c(x[[nm]],y[[nm]])
    }
  }
  com.rdl <- RangesList(.outlist)
  c(com.rdl, x[x.only.names], y[y.only.names])
}
