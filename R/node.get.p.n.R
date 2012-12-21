node.get.p.n <-
function(anode, win.size) {
  res <- list(n=NA, p=NA)
  if(!is.leaf(anode)) {
    ewins <- as.numeric(labels(anode))/win.size
    n <- (max(ewins)-min(ewins))
    r <- attr(anode, 'members')
    p <- r/n
    res <- list(n=n, p=p)
  }
  res
}
