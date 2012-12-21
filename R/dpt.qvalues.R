dpt.qvalues <-
function(flat.test.report,
                        cutoff.by.winsize=T,
                        winsize=50,
                        cutoff.by.winsize.pct=0.9,
                        pvalue.var="p.value.1",
                        split.vars=c("pattern", "contrast")) {
  
### handle NA pvalues
### flat.test.report is a RangedData class
  qvalues.init <- rep(NA, nrow(flat.test.report))
  pval.na.idx <- is.na(flat.test.report[[pvalue.var]])  

### split by pattern, contrast, and possibly winsize
  cat("# of sites:", nrow(flat.test.report), "\n")
  cat("# of NA pvalues:", sum(pval.na.idx), "\n")
  
  in.sub <- flat.test.report[!pval.na.idx,]
  pvals <- in.sub[[pvalue.var]]
  in.width <- width(in.sub)
  
  fac.list <- list(patt=in.sub$pattern, ctr=in.sub$contrast)

  pvals.split <- split(data.frame(p=pvals, w=in.width), fac.list)
  qvals.split <- lapply(pvals.split,
                        function(x, cut.pct) {
                         
                          if(nrow(x)==0) {
                            return(numeric())
                          } else {
                            w <- x$w
                            p <- x$p
                            w.fac.max <- quantile(w, cut.pct)
                            w[w>w.fac.max] <- w.fac.max+1
                            w.fac <- factor(w)
                            p.split <- split(p, w.fac)
                            unsplit(lapply(p.split, compute.qvalues, lambda=0), w.fac)
                          }
                        },
                        cut.pct=cutoff.by.winsize.pct)
  
### unsplit the calculated qvalues
  qvalues <- unsplit(qvals.split, fac.list)
  qvalues.init[!pval.na.idx] <- qvalues
  qvalues.init
}
