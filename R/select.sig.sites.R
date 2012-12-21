#' @export
select.sig.sites <-
function(flat.test.report,
                             method.p.sel="qvalue",
                             cutoff.p.sel=0.05,
                             cutoff.v.sel=NA,
                             winsize=50,
                             cutoff.by.winsize=T,
                             cutoff.by.winsize.pct=0.9,
                             pvalue.var="p.value.1",
                             value.var="value.1",
                             split.vars=c("pattern", "contrast"),
                             cols.remove=c("space","seq", "ID")) {
  
  method.p.sel <- match.arg(method.p.sel)
  if(method.p.sel == "qvalue") {
    .rd.full <- trans2rangedData(dpt.report.postprocess(flat.test.report,
                                                        cutoff.by.winsize=cutoff.by.winsize,
                                                        winsize=winsize,
                                                        cutoff.by.winsize.pct=cutoff.by.winsize.pct,
                                                        pvalue.var="p.value.1",
                                                        ##value.var="value.1",
                                                        split.vars=split.vars,
                                                        cols.remove=cols.remove),
                                 ext.end=NULL)
    
    .rd.sel <- cutoff.sites(.rd.full,
                            qvalue.cutoff=cutoff.p.sel,
                            value.cutoff=cutoff.v.sel,
                            value.var=value.var,
                            qvalue.var="qvalue")
    .df.sel <- as.data.frame(.rd.sel)
    .df.sel <- .df.sel[,!(names(.df.sel)%in%c("space"))]
    
  } else {stop("Cutoff method not implemented!")}
  
  attr(.df.sel, 'testing') <- list(method=method.p.sel, cutoff=cutoff.p.sel)
  .df.sel
}
