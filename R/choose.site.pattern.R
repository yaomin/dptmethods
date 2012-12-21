choose.site.pattern <-
function(e.res, start, end, e.cutoff=-1, site.cutoff=.66) {
  ##browser()
  site.patt.df <- subset(e.res,
                         subset=as.integer(row.names(e.res))>=as.integer(start) & as.integer(row.names(e.res)) <= as.integer(end)) > e.cutoff
  patt.counts <- colSums(site.patt.df)
  n.patt <- nrow(site.patt.df)
  names.patt <- dimnames(site.patt.df)[[2]]
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
                   best.valid=patt.best.valid)
  out.list
}
