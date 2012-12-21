sel.sigsites <-
function(flat.test.report,
                         pattern,
                         contrast,
                         method =c(
                           "holm",
                           "hochberg",
                           "hommel",
                           "bonferroni",
                           "BH",
                           "BY",
                           "fdr",
                           "none",
                           "qvalue"),
                         pvalue.var="p.value.1",
                         cutoff=0.05)
{
  ## See Table below from flat.test.report on columns pattern and contrast
  ## table(test.report.flat[,c('pattern','contrast')])
  ##          contrast
  ## pattern 011_0-11 001_-110 010_-101 010_001
  ##     011     5780        0        0       0
  ##     001        0      749        0       0
  ##     100     1783        0        0       0
  ##     110        0     4639        0       0
  ##     010        0        0      827       0
  ##     101        0        0     4112       0
  ##     111    30992    30992    30992   30992
  ##     **1        0        0        0   36847
  ##     *1*        0        0        0   36719
  ##     1**        0        0        0   36947

  ##browser()
  method <- match.arg(method)
  
  idx.patt <- flat.test.report$pattern %in% pattern
  idx.contr <- flat.test.report$contrast %in% contrast
  idx.start <- idx.patt & idx.contr

  idx.out <- NULL
  if(sum(idx.start)>0) {
    pvals <- subset(flat.test.report,
                    subset=idx.start,
                    select=pvalue.var,
                    drop=T)
    if(method != "qvalue") {
      pvals.adj <- p.adjust(pvals, method=method)
    } else if(method == "qvalue") {
      pvals.adj <- compute.qvalues(pvals,lambda=0)
    } else stop("Not supported multiple testing method provided!")
    
    idx.out <- which(idx.start)[pvals.adj < cutoff]
  }
  idx.out
}
