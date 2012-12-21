select.sites <-
function(reports,
                         cut.fdr=0.05, cut.effect = 1,
                         cut.base=0.2, cut.eq=0.2,
                         col.select=NULL, win.size=100) {

  if(is.null(col.select)) col.select <- seq(ncol(reports))
  sel.fdr <- p.adjust(reports$pvalue, 'fdr') < cut.fdr
  sel.effect <- abs(reports$effect) > cut.effect
  sel.base <- reports$pbase > cut.base
  sel.eq <- reports$peq > cut.eq

  sel <- apply(cbind(sel.fdr, sel.effect, sel.base, sel.eq), 1, all)
  res <- reports[sel,col.select]
  res$end <- res$end+win.size
  size <- res$end-res$start
  res <- cbind(res[,c(1,2)],size,res[,-c(1,2)])
  reorder <- order(res$pvalue)
  res[reorder,]
}
