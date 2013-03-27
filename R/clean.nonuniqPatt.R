#' @export
clean.nonuniqPatt <-
function(sitelist) {
  l.n <- length(sitelist)
  idx <- common.expandIdx(l.n, l.n, self.included=F)
  l.nms <- names(sitelist)
  for(i in seq(nrow(idx))) {
    nr.before <- unlist(lapply(sitelist[idx[i,]], nrow))
    sitelist[idx[i,]] <- compete.patt(sitelist[idx[i,]])
    nr.after <- unlist(lapply(sitelist[idx[i,]], nrow))
    cat("comparing",l.nms[idx[i,]],":",
        nr.before,",",nr.after,
        "\t;", nr.after-nr.before,
        "\n")
  }
  sitelist
}
