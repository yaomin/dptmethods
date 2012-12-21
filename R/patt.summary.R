patt.summary <-
function(a.values, site.cutoff=0.66) {
  ##browser()
  ##a.values <- as.data.frame(values(a.rd))[,-1]
  patt.counts <- colSums(a.values)
  patt.counts.ord.12 <- sort(patt.counts, decreasing=T)[1:2]
  patt.diff.score <- abs(diff(patt.counts.ord.12))*2/sum(patt.counts.ord.12)
  n.patt <- nrow(a.values)
  names.patt <- dimnames(a.values)[[2]]
  patts.exist <- names.patt[patt.counts>0]
  patt.best <- which.max(patt.counts)
  patt.best.name <- names.patt[patt.best]
  patt.cutoff <- floor(n.patt*site.cutoff)
  patt.best.valid <- patt.counts[patt.best] > patt.cutoff
  out.list <- list(counts=patt.counts,
                   N=n.patt,
                   patterns=paste(patts.exist, collapse=","),
                   best.pattern=patt.best.name,
                   best.cutoff=patt.cutoff,
                   best.valid=patt.best.valid,
                   diff.score=patt.diff.score)
  out.list
}
