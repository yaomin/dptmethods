select.sig.sites.2 <-
function(flat.test.report, contrast.objs, method.sel="fdr", cutoff.p.sel=0.05,
                               winsize=50, cutoff.by.winsize=T, ul.by.winsize=30)
{

  ### remove NA pvalues
  pval.na.idx <- which(is.na(flat.test.report$p.value.1))
  if(length(pval.na.idx) > 0) {
    cat("# of missing p values: ", length(pval.na.idx),"\n")
    flat.test.report <- flat.test.report[-pval.na.idx,]
  }
  ### process selection
  
  out <- NULL
  if(cutoff.by.winsize) {
    wsizes <- as.integer(names(with(flat.test.report, table(end-start))))
    wsize.ul <- ul.by.winsize * winsize
    singles <- wsizes[wsizes<=wsize.ul]
    for(asize in singles) {
      a.report <- subset(flat.test.report, subset=(end-start)==asize)
      idx.sel <- NULL
      for(i in seq(along=contrast.objs)) {
        idx.sel <- union(idx.sel, sel.sigsites(a.report,
                                               pattern = contrast.objs[[i]]$pattern,
                                               contrast = contrast.objs[[i]]$contrast,
                                               method = method.sel,
                                               pvalue.var = contrast.objs[[i]]$pvalue.var,
                                               cutoff = cutoff.p.sel)
                         )
      }
      out <- rbind(out, a.report[idx.sel,])
    }
    a.report <- subset(flat.test.report, subset=(end-start)>wsize.ul)
    idx.sel <- NULL
    for(i in seq(along=contrast.objs)) {
      idx.sel <- union(idx.sel, sel.sigsites(a.report,
                                             pattern = contrast.objs[[i]]$pattern,
                                             contrast = contrast.objs[[i]]$contrast,
                                             method = method.sel,
                                             pvalue.var = contrast.objs[[i]]$pvalue.var,
                                             cutoff = cutoff.p.sel)
                       )
    }
    out <- rbind(out, a.report[idx.sel,])
  } else {
    idx.sel <- NULL
    for(i in seq(along=contrast.objs)) {
      idx.sel <- union(idx.sel, sel.sigsites(flat.test.report,
                                             pattern = contrast.objs[[i]]$pattern,
                                             contrast = contrast.objs[[i]]$contrast,
                                             method = method.sel,
                                             pvalue.var = contrast.objs[[i]]$pvalue.var,
                                             cutoff = cutoff.p.sel)
                       )
    }
    out <- test.report.flat[idx.sel,]
  }
  attr(out, 'contrasts') <- contrast.objs
  attr(out, 'testing') <- list(method=method.sel, cutoff=cutoff.p.sel)
  out
}
