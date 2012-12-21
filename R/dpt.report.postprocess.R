dpt.report.postprocess <-
function(flat.test.report,
                                   cutoff.by.winsize=T,
                                   winsize=50,
                                   cutoff.by.winsize.pct=0.9,
                                   pvalue.var="p.value.1",
                                   split.vars=c("pattern", "contrast"),
                                   cols.remove=c("space","seq", "ID")) {
  ##browser()
  .dt <- trans2rangedData(remove.dupsite(flat.test.report), ext.end=winsize-1)
  .dt$qvalue <- dpt.qvalues(.dt,
                            cutoff.by.winsize=cutoff.by.winsize,
                            winsize=winsize,
                            cutoff.by.winsize.pct=cutoff.by.winsize.pct,
                            pvalue.var=pvalue.var,
                            split.vars=split.vars)
  .dt <- remove.dupcontr(.dt, pvalue.var=pvalue.var)
  .dt[,!(names(.dt)%in%cols.remove)]
}
