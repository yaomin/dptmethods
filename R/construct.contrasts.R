#'  @export
construct.contrasts <-
function(patts) {
  ##browser()
  patt.ncol <- ncol(as.matrix(patts))
  patt.sum <- colSums(patts)
  out <- vector("list", patt.ncol)
  for(i in seq(patt.ncol)) {
    ipatt <- patts[,i]
    if(patt.sum[i] != nrow(patts)) out[[i]] <- default.contrastList(ipatt)
    else {
      out[[i]] <- default.contrastList(patts[,patt.sum!=nrow(patts)])
    }
  }
  out
}
