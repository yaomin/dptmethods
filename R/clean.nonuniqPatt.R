clean.nonuniqPatt <-
function(sitelist) {
  ##browser()
  l.n <- length(sitelist)
  idx <- common.expandIdx(l.n, l.n, self.included=F)
  l.nms <- names(sitelist)
  for(i in seq(nrow(idx))) {
    sitelist[idx[i,]] <- compete.patt(sitelist[idx[i,]])
    cat("comparing",l.nms[idx[i,]],"\n")
  }
  sitelist
}
